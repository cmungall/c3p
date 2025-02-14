"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is defined as a naturally occurring, basic nitrogen compound (mostly heterocyclic).
    Amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars, and antibiotics are
    not normally regarded as alkaloids. Compounds where nitrogen is exocyclic (e.g., dopamine, mescaline,
    serotonin) are usually classed as amines rather than alkaloids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not n_atoms:
        return False, "No nitrogen atoms found"

    # Check for basic nitrogen atoms
    basic_nitrogens = []
    for atom in n_atoms:
        if atom.GetFormalCharge() != 0:
            continue  # Exclude charged nitrogens (quaternary ammonium)
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            # SP3 hybridized nitrogen is often basic
            basic_nitrogens.append(atom)
        elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            # For SP2 hybridized nitrogen, check if it's pyridine-like
            num_h = atom.GetTotalNumHs()
            if num_h == 1:
                basic_nitrogens.append(atom)

    if not basic_nitrogens:
        return False, "No basic nitrogen atoms found"

    # Exclude amino acids and peptides by detecting peptide bonds
    peptide_bond = Chem.MolFromSmarts("N[C;D2](=O)C")  # N-C(=O)-C pattern
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Contains peptide bond, may be a peptide or protein"

    # Exclude nucleotides and nucleic acids by detecting nucleobases
    nucleic_acid_bases = [
        Chem.MolFromSmarts("c1nc[nH]c(=O)[nH]1"),  # Cytosine
        Chem.MolFromSmarts("c1cc(=O)[nH]c(=O)[nH]1"),  # Uracil/Thymine
        Chem.MolFromSmarts("c1ncnc2ncnn12"),  # Adenine
        Chem.MolFromSmarts("c1[nH]c2c(n1)nc(nc2)N"),  # Guanine
    ]
    for base in nucleic_acid_bases:
        if mol.HasSubstructMatch(base):
            return False, "Contains nucleic acid base"

    # Exclude amino sugars by detecting sugar rings with amino groups
    amino_sugar = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](N)[C@@H](O)[C@H](O)[C@H]1O")  # Simple amino sugar pattern
    if mol.HasSubstructMatch(amino_sugar):
        return False, "Contains amino sugar moiety"

    # Exclude antibiotics - often contain characteristic lactone or glycopeptide structures
    antibiotic_patterns = [
        Chem.MolFromSmarts("C(=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # General lactone ring
        Chem.MolFromSmarts("C1=CNC(=O)N=C1"),  # Beta-lactam ring
    ]
    for pattern in antibiotic_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains antibiotic-like structure"

    # Exclude certain functional groups associated with excluded classes
    excluded_groups = [
        Chem.MolFromSmarts("P(=O)(O)O"),  # Phosphate group in nucleotides
        Chem.MolFromSmarts("C=O"),        # Carbonyl group in amino acids and peptides
        Chem.MolFromSmarts("C#N"),        # Nitrile group
    ]
    for group in excluded_groups:
        if mol.HasSubstructMatch(group):
            group_name = group.GetProp("_Name") if mol.HasProp("_Name") else "Excluded group"
            return False, f"Contains {group_name}, may be an excluded compound"

    # Check if nitrogen is part of an amine group (exocyclic)
    amine_n = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary or secondary amine nitrogen
    if mol.HasSubstructMatch(amine_n):
        # Exclude if nitrogen is exocyclic and not in a ring
        for match in mol.GetSubstructMatches(amine_n):
            atom = mol.GetAtomWithIdx(match[0])
            if not atom.IsInRing():
                return False, "Nitrogen is exocyclic, molecule may be an amine rather than an alkaloid"

    # At this point, we have a molecule with at least one basic nitrogen atom
    # that is not excluded by the other criteria
    return True, "Molecule is classified as an alkaloid"
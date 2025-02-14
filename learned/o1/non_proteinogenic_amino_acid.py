"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:64049 non-proteinogenic amino acid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is any amino acid that is not among the standard amino acids
    encoded by the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens (necessary for accurate matching)
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for amine and carboxylic acid groups
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary or secondary amine
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")  # Carboxylic acid group

    # Check for amine group
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No primary or secondary amine group found"

    # Check for carboxylic acid group
    acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # Exclude peptides by checking for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C;!$([C;!H0])]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if peptide_bonds:
        return False, "Peptide bonds detected, molecule may be a peptide"

    # Generate canonical SMILES for the molecule
    mol_canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    # List of standard amino acids canonical SMILES (with stereochemistry)
    standard_amino_acids_smiles = [
        "N[C@@H](C)C(=O)O",                          # L-Alanine
        "N[C@@H](CC(=O)O)C(=O)O",                    # L-Aspartic acid
        "N[C@@H](C(=O)N)C(=O)O",                     # L-Asparagine
        "N[C@@H](CS)C(=O)O",                         # L-Cysteine
        "NCC(=O)O",                                  # Glycine
        "N[C@@H](CCC(=O)O)C(=O)O",                   # L-Glutamic acid
        "N[C@@H](CC(=O)N)C(=O)O",                    # L-Glutamine
        "N[C@@H](Cc1c[nH]cn1)C(=O)O",                # L-Histidine
        "N[C@@H](I)C(=O)O",                          # L-Isoleucine
        "N1[C@@H]([C@H](CCC1)C(=O)O)C(=O)O",         # L-Proline
        "N[C@@H](CC(C)C)C(=O)O",                     # L-Leucine
        "N[C@@H](CO)C(=O)O",                         # L-Serine
        "N[C@@H](Cc1ccccc1)C(=O)O",                  # L-Phenylalanine
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",               # L-Tyrosine
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",          # L-Tryptophan
        "N[C@@H](CCCN)C(=O)O",                       # L-Lysine
        "N[C@@H](CCS)C(=O)O",                        # L-Methionine
        "N[C@@H](C(C)O)C(=O)O",                      # L-Threonine
        "N[C@@H](C(C)C)C(=O)O",                      # L-Valine
        "N[C@@H](CSO)C(=O)O",                        # L-Selenocysteine
        "N[C@@H](C(O)O)C(=O)O",                      # L-Serine
        # Include D-forms if necessary
    ]

    # Generate canonical SMILES for standard amino acids
    standard_amino_acids_canon_smiles = set()
    for smi in standard_amino_acids_smiles:
        std_mol = Chem.MolFromSmiles(smi)
        if std_mol:
            std_mol = Chem.AddHs(std_mol)
            std_mol_canon = Chem.MolToSmiles(std_mol, isomericSmiles=True)
            standard_amino_acids_canon_smiles.add(std_mol_canon)

    # Check if the molecule matches any standard amino acid
    if mol_canonical_smiles in standard_amino_acids_canon_smiles:
        return False, "Molecule is a standard (proteinogenic) amino acid"

    # Alternatively, compare InChIKeys
    mol_inchi_key = Chem.MolToInchiKey(mol)
    for smi in standard_amino_acids_smiles:
        std_mol = Chem.MolFromSmiles(smi)
        if std_mol:
            std_mol = Chem.AddHs(std_mol)
            std_inchi_key = Chem.MolToInchiKey(std_mol)
            if mol_inchi_key == std_inchi_key:
                return False, "Molecule is a standard (proteinogenic) amino acid"

    return True, "Molecule is a non-proteinogenic amino acid (contains amine and carboxylic acid groups but is not a standard amino acid)"
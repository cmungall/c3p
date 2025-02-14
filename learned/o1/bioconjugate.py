"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is a molecular entity consisting of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for biological molecules
    amino_acid_pattern = Chem.MolFromSmarts("N[C;!H0][CX3](=O)[O,N]")
    sugar_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1O")  # Simplified glucose pattern
    nucleotide_pattern = Chem.MolFromSmarts("n1c([nH])cnc1")  # Purine/pyrimidine ring
    lipid_pattern = Chem.MolFromSmarts("C(=O)O[C;H2][C;H2]")  # Ester-linked fatty acid
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;!H0]=O")  # Peptide bond

    # Find matches
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    nucleotide_matches = mol.GetSubstructMatches(nucleotide_pattern)
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Initialize count
    bio_units = 0

    # Use a set to store the indices of matched atoms to avoid double-counting
    matched_atoms = set()

    # Check for amino acids
    for match in amino_acid_matches:
        if not matched_atoms.intersection(match):
            bio_units += 1
            matched_atoms.update(match)

    # Check for sugars
    for match in sugar_matches:
        if not matched_atoms.intersection(match):
            bio_units += 1
            matched_atoms.update(match)

    # Check for nucleotides
    for match in nucleotide_matches:
        if not matched_atoms.intersection(match):
            bio_units += 1
            matched_atoms.update(match)

    # Check for lipids
    for match in lipid_matches:
        if not matched_atoms.intersection(match):
            bio_units += 1
            matched_atoms.update(match)

    # Optionally, check for peptide bonds as connectors
    has_peptide_bond = len(peptide_matches) > 0

    # Determine if the molecule is a bioconjugate
    if bio_units >= 2:
        return True, f"Contains {bio_units} biological units covalently linked"
    else:
        return False, "Does not contain at least 2 covalently linked biological units"
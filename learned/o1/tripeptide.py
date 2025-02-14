"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:16699 tripeptide
"""

from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide consisting of three amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the peptide bond pattern (amide bond): O=C-N
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    if num_peptide_bonds != 2:
        return False, f"Found {num_peptide_bonds} peptide bonds, need exactly 2 for tripeptide"

    # Define the amino acid residue pattern in a peptide chain
    # Pattern: N-C-C(=O)
    residue_pattern = Chem.MolFromSmarts("[N][CH1][C](=O)")
    residue_matches = mol.GetSubstructMatches(residue_pattern)
    num_residues = len(residue_matches)

    if num_residues != 3:
        return False, f"Found {num_residues} amino acid residues, need exactly 3 for tripeptide"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"
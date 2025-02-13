"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide that consists of three amino acid residues connected by peptide linkages.

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

    # Pattern for alpha amino acids - more generalized to include non-chiral carbons
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4](C)C=O")
    if alpha_carbon_pattern is None:
        return False, "Failed to create alpha carbon pattern"

    # Peptide bond pattern: -C(=O)N-, generalized more to capture broader peptide linkage configuration
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if peptide_bond_pattern is None:
        return False, "Failed to create peptide bond pattern"

    # Finding matches
    alpha_matches = len(mol.GetSubstructMatches(alpha_carbon_pattern))
    peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))

    # Tripeptides have exactly 3 residues which translates to 2 peptide bonds
    if alpha_matches >= 3 and peptide_bonds == 2:
        return True, "Contains three alpha carbons with two peptide linkages indicating a tripeptide"
    else:
        return False, f"Contains {alpha_matches} alpha carbons and {peptide_bonds} peptide bonds, insufficient or excess for a tripeptide"
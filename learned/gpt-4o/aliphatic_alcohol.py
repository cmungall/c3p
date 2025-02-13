"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of hydroxyl group attachment (C-OH)
    oh_pattern = Chem.MolFromSmarts("[CX4H2,CX4H,CX4]-[OH]")  # Saturated carbon with hydroxyl
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No aliphatic hydroxyl group found"

    # Verify the absence of aromaticity
    if mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return False, "Contains aromatic rings; not purely aliphatic"

    return True, "Contains aliphatic hydroxyl group(s) on a non-aromatic carbon chain"
"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxy groups connected by a carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dicarboxylic acid, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for dicarboxylic acid pattern
    dicarboxylic_pattern = Chem.MolFromSmarts("[C](=O)(O)[C]([#6])~[#6][C](=O)(O)")
    dicarboxylic_matches = mol.GetSubstructMatches(dicarboxylic_pattern)

    # Exclude molecules with more than two carboxyl groups
    carboxyl_pattern = Chem.MolFromSmarts("[C](=O)(O)")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 2:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, expected 2"

    if dicarboxylic_matches:
        return True, "Contains two carboxyl groups connected by a carbon chain"
    else:
        return False, "Does not match the dicarboxylic acid pattern"
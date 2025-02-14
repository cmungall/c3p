"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains exactly two carboxylic acid (COOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Get matches for the carboxylic acid pattern
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Check the number of carboxylic acid groups
    if len(carboxylic_acid_matches) == 2:
        return True, "Contains exactly two carboxylic acid groups"
    else:
        return False, f"Contains {len(carboxylic_acid_matches)} carboxylic acid groups, not exactly 2"
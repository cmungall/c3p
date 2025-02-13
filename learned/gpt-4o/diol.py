"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol contains exactly two hydroxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find hydroxyl groups (–OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check number of hydroxyl groups
    if len(hydroxy_matches) != 2:
        return False, f"Found {len(hydroxy_matches)} hydroxy groups, need exactly 2"
    
    # The molecule has two hydroxy groups, classify as diol
    return True, "Contains exactly two hydroxy groups"

# Example usage:
# result, reason = is_diol("OC(CO)CO")  # Should return True, "Contains exactly two hydroxy groups"
# print(result, reason)
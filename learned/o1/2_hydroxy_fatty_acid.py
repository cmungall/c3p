"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:136568 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carboxylic acid group
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H]')
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)

    if not carboxy_matches:
        return False, "No carboxylic acid group found"
    
    # Define SMARTS pattern for alpha-hydroxy group adjacent to carboxylic acid
    alpha_hydroxy_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H]-[CX4]-[OX2H]')
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    
    if alpha_hydroxy_matches:
        return True, "Contains carboxylic acid group with hydroxy on alpha carbon (2-hydroxy fatty acid)"
    else:
        return False, "No hydroxy group found on alpha carbon adjacent to carboxylic acid group"
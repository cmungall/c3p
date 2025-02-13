"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of dimethylisoalloxazine with a substituent on the 10th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    # Convert SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for dimethylisoalloxazine core
    pteridine_pattern = Chem.MolFromSmarts("c1cc2ncnc-3[nH]c(=O)n(c3=O)c2cc1")
    
    # Check pteridine core
    if not mol.HasSubstructMatch(pteridine_pattern):
        return False, "Core pteridine structure missing"
    
    # Define SMARTS pattern to identify dimethyl groups at specific positions
    dimethyl_pattern = Chem.MolFromSmarts("c1c([CH3])c([CH3])cc1")
    
    # Check for dimethyl positions
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Required dimethyl groups missing"
    
    # Define a pattern for substitution at the 10th position
    ten_position_pattern = Chem.MolFromSmarts("n1cnc2c1nc(=O)[nH]c2=O")
    
    # Check substitution at the 10th position
    if not mol.HasSubstructMatch(ten_position_pattern):
        return False, "No substitution at the 10th position"
    
    return True, "Valid flavin structure with dimethylisoalloxazine core and 10th position substitution"
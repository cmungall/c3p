"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is defined as a sesquiterpene lactone with a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for decalin or cyclodecane-like structure (10-membered carbon-containing ring)
    cyclodecane_pattern = Chem.MolFromSmarts("C1CCCCCCCC1")
    if not mol.HasSubstructMatch(cyclodecane_pattern):
        return False, "No cyclodecane-like structure found"

    # Check for lactone group (-C(=O)O-)
    lactone_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"

    # Optionally, check for decoration like hydroxyl groups, double bonds
    # This can be expanded based on the specific chemotype details

    return True, "Contains cyclodecane skeleton with lactone group typical of germacranolides"
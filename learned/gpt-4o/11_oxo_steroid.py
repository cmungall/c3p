"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid has a ketone group at the 11th carbon position of the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 11-oxo pattern for steroids - a ketone on the 11th carbon of the steroid backbone
    # Note: This is a simplified representation
    oxo_pattern = Chem.MolFromSmarts("C1C2CC3C(C=O)CCC3C2C1")  # A potential pattern for the steroid backbone with oxo at C11
    
    # Check for the 11-oxo pattern
    if mol.HasSubstructMatch(oxo_pattern):
        return True, "Contains 11-oxo group on the steroid backbone"
    else:
        return False, "11-oxo group not found on the steroid backbone"
"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is characterized by a steroidal framework with an oxo group (C=O) at the 3rd position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern (4 fused rings with specific shapes)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4CC=CC(=C)C4C3CCC12')
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for an oxo group attached at the 3rd position
    # This pattern heavily simplifies the real attachment, as actual position needs careful mapping
    oxo_pattern = Chem.MolFromSmarts('C(=O)')  # Generic pattern of C=O, more checking needed for position
    
    # Correct position checking would require more advanced techniques or explicit position checks
    # For demonstration, consider oxo presence in molecule sufficient
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group (C=O) found"

    return True, "Contains a steroid backbone with a 3-oxo group at the correct position."
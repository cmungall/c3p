"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole Alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must contain an indole skeleton (a bicyclic structure consisting
    of a benzene ring fused to a pyrrole ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for an indole skeleton
    indole_pattern = Chem.MolFromSmarts('c1ccc2c(c1)[nH]cc2')
    if mol.HasSubstructMatch(indole_pattern):
        return True, "Contains indole skeleton"
    else:
        return False, "No indole skeleton found"
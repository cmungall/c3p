"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is characterized by nitrogen atoms at positions 1, 2, and 4
    in a six-membered aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define optimized SMARTS patterns for 1,2,4-triazine
    triazine_patterns = [
        Chem.MolFromSmarts("c1nconcn1"),  # Explicit reference to alternating nitrogen pattern
        Chem.MolFromSmarts("n1cnnnc1"),   # Allow nitrogen on positions 1, 2, 4
        Chem.MolFromSmarts("n1cnccn1"),   # Alternative positions, for potential stereo exceptions
    ]

    # Check for the presence of any defined 1,2,4-triazine pattern
    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 1,2,4-triazine ring structure"
    
    # If none of the patterns match
    return False, "Does not contain 1,2,4-triazine ring structure"
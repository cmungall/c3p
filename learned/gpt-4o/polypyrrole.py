"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is defined as a compound with two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pyrrole pattern (5-membered ring with NH)
    pyrrole_pattern = Chem.MolFromSmarts("n1cccc1")

    # Find matches for pyrrole units
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Check if there are two or more pyrrole units
    if len(pyrrole_matches) >= 2:
        return True, f"Contains {len(pyrrole_matches)} pyrrole units"
    
    return False, "Less than 2 pyrrole units found"
"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a three-membered ring containing an oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for epoxide (three-membered ring with one oxygen)
    epoxide_pattern = Chem.MolFromSmarts("[C]1[C][O]1")

    # Check for a match of the pattern
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains a three-membered ring with an oxygen"

    return False, "Does not contain an epoxide ring"
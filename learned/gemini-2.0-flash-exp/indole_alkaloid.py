"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton.

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

    # Look for indole substructure
    indole_pattern = Chem.MolFromSmarts("c1c[nH]c2ccccc12")
    if mol.HasSubstructMatch(indole_pattern):
        return True, "Contains an indole substructure"
    else:
        return False, "Does not contain an indole substructure"
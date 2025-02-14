"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is an aldehyde where the carbonyl group is directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring carbon directly attached to a carbonyl carbon
    # The aromatic carbon is [c]
    # The carbonyl carbon is [CX3]=O
    # The connection is a single bond
    aldehyde_pattern = Chem.MolFromSmarts("[c][CX3](=[OX1,O])")


    if not mol.HasSubstructMatch(aldehyde_pattern):
       return False, "No aldehyde group directly attached to aromatic ring found"

    return True, "Contains an aromatic ring with a directly attached aldehyde group."
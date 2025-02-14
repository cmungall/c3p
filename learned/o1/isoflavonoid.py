"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is any molecule containing the isoflavone core structure,
    which is a 3-phenylchromen-4-one.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoflavone core SMARTS pattern
    # The isoflavone core: 3-phenylchromen-4-one
    isoflavone_smarts = 'O=C1C=CC2=CC=CC=C2O1-c3ccccc3'

    isoflavone_pattern = Chem.MolFromSmarts(isoflavone_smarts)
    if isoflavone_pattern is None:
        return False, "Invalid SMARTS pattern for isoflavone core"

    # Search for the isoflavone core
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Contains isoflavone core structure"
    else:
        return False, "Does not contain isoflavone core structure"
"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is defined as a compound containing a 2-aryl-1-benzopyran-4-one
    (2-arylchromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the flavone core
    # Flavone core: chromone (1-benzopyran-4-one) with an aryl group at position 2
    flavone_smarts = 'O=C1C=CC(=C1C2=CC=CC=C2)'

    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    if flavone_pattern is None:
        return False, "Invalid SMARTS pattern for flavone core"

    # Check if the molecule matches the flavone pattern
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains flavone core structure (2-arylchromen-4-one skeleton)"
    else:
        return False, "Does not contain flavone core structure"
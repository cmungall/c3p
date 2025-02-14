"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones

Determines if a molecule is a flavanone based on its SMILES string.
A flavanone has a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavanone core SMARTS pattern
    # Flavanone core: 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one
    flavanone_pattern = Chem.MolFromSmarts('O=C1CCOc2ccccc12')

    if flavanone_pattern is None:
        return None, "Error in defining flavanone SMARTS pattern"

    # Check for flavanone core
    if mol.HasSubstructMatch(flavanone_pattern):
        return True, "Contains flavanone core structure"
    else:
        return False, "Does not contain flavanone core structure"
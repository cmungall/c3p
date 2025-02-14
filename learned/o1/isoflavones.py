"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton (3-aryl-4H-chromen-4-one)
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Corrected SMILES for the isoflavone core
    isoflavone_smiles = 'O=C1c2ccccc2Oc3ccccc13'  # Isoflavone core structure
    isoflavone_core = Chem.MolFromSmiles(isoflavone_smiles)
    if isoflavone_core is None:
        return False, "Error in isoflavone core SMILES"

    # Check if the molecule contains the isoflavone core
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "Does not contain the 3-aryl-1-benzopyran-4-one skeleton"

    # The molecule contains the isoflavone core
    return True, "Contains the 3-aryl-1-benzopyran-4-one skeleton characteristic of isoflavones"
"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    This classification requires the presence of a hydroxyl group at the 17-alpha position
    in conjunction with the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid backbone SMARTS pattern
    steroid_backbone_smarts = "C1=CC2CC3CCC4CCCC(C)(C)O3CC2C(C1)=O"

    # Define SMARTS pattern for 17alpha-hydroxy group in a simplified way
    # The 'alpha' hydroxy group could be represented in a way specific to steroid hydroxylation patterns
    # Here, we use a general pattern for hydroxy over a four-ring structure, though it's a simplistic representation
    hydroxy_17alpha_smarts = "[C@H]1([C@@H](C2CCCC3CC(O)CCCC3=O)O)CCC2C1"

    # Check for the steroid backbone
    steroid_backbone = Chem.MolFromSmarts(steroid_backbone_smarts)
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for presence of 17alpha-hydroxy group
    hydroxy_17alpha = Chem.MolFromSmarts(hydroxy_17alpha_smarts)
    if not mol.HasSubstructMatch(hydroxy_17alpha):
        return False, "No 17alpha-hydroxy group found"

    return True, "Contains steroid backbone with 17alpha-hydroxy group"
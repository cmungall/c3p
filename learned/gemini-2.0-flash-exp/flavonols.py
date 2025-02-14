"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a flavone with a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core pattern using SMARTS.
    # The core structure: C1=CC=C(C=C1)-C2=C(C(=O)C3=CC=CC=C3O2)-O is a good starting point.
    # But we will be more specific with query.
    # It's crucial to define the aromaticity of the rings in this pattern and include the O in position 1 (or 4) of the pyran ring.
    # We also specify the correct atom connectivities
    flavone_core = Chem.MolFromSmarts('c1ccccc1-c2oc(c3ccccc3)c(=O)c2-O') # 1-2-3
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Not a flavone core structure"

    # Check for the hydroxyl group at position 3 of the pyran ring (C ring).
    # This is the carbon directly attached to the carbonyl group in the flavone core
    # We will try to identify the C next to the carbonyl C, and check it has -OH group.
    hydroxyl_at_3 = Chem.MolFromSmarts('c1ccccc1-c2oc(-[C;H1]([O;H1])c(=O)c3ccccc3)c2') #1-2-3-OH
    if not mol.HasSubstructMatch(hydroxyl_at_3):
      return False, "No hydroxyl group at position 3 of the flavone core"

    # If both the flavone core and 3-OH are found, it's a flavonol
    return True, "Flavonol: Flavone core with hydroxyl at position 3"
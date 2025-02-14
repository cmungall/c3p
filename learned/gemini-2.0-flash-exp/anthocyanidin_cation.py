"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an aglycon of anthocyanin cation, characterized by a flavylium core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavylium core SMARTS pattern
    flavylium_core_pattern = Chem.MolFromSmarts("[O+]=C1C=CC=C(C2=CC=CC=C2)C1")

    # Check for the flavylium core
    if not mol.HasSubstructMatch(flavylium_core_pattern):
      return False, "Flavylium core not found"

    # Check for -OH groups attached to the core - this is a weak criteria
    oh_group_pattern = Chem.MolFromSmarts("[O;H]-c1ccccc1")
    oh_matches = mol.GetSubstructMatches(oh_group_pattern)

    # Check for sugars attached directly at core
    sugar_pattern_1 = Chem.MolFromSmarts("C[C,O]1[C,O][C,O][C,O][C,O][C,O]1~[OX2]")
    sugar_pattern_2 = Chem.MolFromSmarts("[C,O]1[C,O][C,O][C,O][C,O][C,O]1~[OX2]")
    
    if mol.HasSubstructMatch(sugar_pattern_1) or mol.HasSubstructMatch(sugar_pattern_2):
      return False, "Sugar moiety present directly at the core"

    return True, "Contains the flavylium cation core (aglycone)"
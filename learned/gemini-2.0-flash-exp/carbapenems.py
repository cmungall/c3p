"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are a class of beta-lactam antibiotics with a specific bicyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core carbapenem structure using SMARTS
    # Matches the bicyclic core with a 4-membered nitrogen containing ring
    # attached to a 5 membered ring
    carbapenem_core_smarts = "[N]12[C@]([H])([C@@H](C1=O)[H])([C@@H]1[C@H]2C=C1)"
    carbapenem_core_pattern = Chem.MolFromSmarts(carbapenem_core_smarts)
    if not mol.HasSubstructMatch(carbapenem_core_pattern):
       return False, "Core carbapenem structure not found"

    # Check for a carboxylic acid group (C(=O)O), which most carbapenems have
    acid_group_smarts = "C(=O)O"
    acid_group_pattern = Chem.MolFromSmarts(acid_group_smarts)

    if not mol.HasSubstructMatch(acid_group_pattern):
       return False, "No carboxylic acid group present"
       
    # Check for the presence of a sulfur substituent at position 3 (general feature)
    sulfur_smarts = "[S]~[C]"
    sulfur_pattern = Chem.MolFromSmarts(sulfur_smarts)
    if not mol.HasSubstructMatch(sulfur_pattern):
       return False, "No sulfur substituent found at the 3 position, this is unusual for carbapenems"

    return True, "Contains the core carbapenem structure and an acid group"
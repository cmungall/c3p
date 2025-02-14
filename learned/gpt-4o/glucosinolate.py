"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates are defined by a central C atom bonded to a glycone group via S,
    a sulfonated oxime group via N, and a side-group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Dynamic SMARTS for glycone group (consider stereochemistry flexibly)
    glycone_pattern = Chem.MolFromSmarts("S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O |flex|")
    if not mol.HasSubstructMatch(glycone_pattern):
        return False, "No thioglucoside group found"

    # Pattern for N-sulfooxy group, where N is doubly bonded to central C
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[CX3](=N\\O[S](=O)(=O)[O-])~*")
    if not mol.HasSubstructMatch(sulfonated_oxime_pattern):
        return False, "No sulfonated oxime group connected"

    # Ensure central carbon connectivity accounting for flexibility
    # Central carbon should have both the S-bonded thioglucoside and N-bonded oxime
    central_c_pattern = Chem.MolFromSmarts("[CX3](S)(=N\\O[S](=O)(=O)[O-])~[*]")
    if not mol.HasSubstructMatch(central_c_pattern):
        return False, "Central carbon with required connectivity not found"

    return True, "Molecule contains the defining features of glucosinolate"
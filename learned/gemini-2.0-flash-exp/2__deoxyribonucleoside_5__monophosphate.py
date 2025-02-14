"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside 5'-monophosphate consists of a 2-deoxyribose sugar, a nucleobase
    attached at the 1' position, and a phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 2'-deoxyribose pattern (5-membered ring with specific substitutions)
    # Using SMARTS to find the ring and check for substitution pattern of H in C2 and C1
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1[CH2][C@H]([CH2])[C@@H](O)[O1]")
    
    if not mol.HasSubstructMatch(deoxyribose_pattern):
         return False, "No 2'-deoxyribose sugar found"

    # Define the phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"


    # Check that the phosphate is attached to C5
    phosphate_connection_pattern = Chem.MolFromSmarts("[CX4][OX2][PX4](=[OX1])([OX2][HX1])[OX2][HX1]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_connection_pattern)

    if not phosphate_matches:
         return False, "Phosphate not at 5' position"

    # Check for nucleobase connected to the deoxyribose C1
    nucleobase_connection_pattern = Chem.MolFromSmarts("[NX3;!H0][CX4]1[CH2][C@H]([CH2])[C@@H](O)[O1]")

    if not mol.HasSubstructMatch(nucleobase_connection_pattern):
        return False, "No nucleobase connected to 1' position"

    return True, "Molecule contains a 2'-deoxyribose with a phosphate on the 5' position and a nucleobase at the 1' position"
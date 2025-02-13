"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for 2'-deoxyribose sugar with no hydroxyl at 2' position and a hydroxy group at 3'
    # Also contains attachments at positions where the base and phosphate would attach
    deoxyribose_pattern = Chem.MolFromSmarts("C1(C(O)C(CO[P](=O)(O)O)O[C@H]1[*])")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar backbone found"

    # Check for a phosphate group at 5' location
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-phosphate group found"

    # Check for a nucleobase (basic structure present in nucleotides)
    # This could be improved to match specific purine/pyrimidine bases
    base_pattern = Chem.MolFromSmarts("[#7]-[#6]=[!#1]")
    if not mol.HasSubstructMatch(base_pattern):
        return False, "No nucleobase found"

    return True, "Contains 2'-deoxyribose with 5'-phosphate and nucleobase"
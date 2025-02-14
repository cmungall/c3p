"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:18115 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.

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

    # Look for 2'-deoxyribo sugar pattern
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]([C@H]([C@@H]([C@H](CO)O)O)O)[C@@H]1C(=O)N")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Look for phosphate group at 5' position
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O[C@@H]1[C@H]([C@@H]([C@H](O1)O)O)N")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-phosphate group found"

    # Check for nucleobase
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # Adenine
        Chem.MolFromSmarts("Nc1ccnc2ncnc12"),  # Guanine
        Chem.MolFromSmarts("Nc1ccnc(=O)n1"),  # Cytosine
        Chem.MolFromSmarts("Oc1ccnc(=O)n1"),  # Uracil
        Chem.MolFromSmarts("Cc1c[nH]c(=O)[nH]c1=O"),  # Thymine
    ]
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No valid nucleobase found"

    return True, "Contains 2'-deoxyribose sugar with 5'-phosphate group and a nucleobase"
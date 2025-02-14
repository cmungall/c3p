"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside contains a pyrimidine base linked to a 2-deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pyrimidine base (generic) with a N at position 1 and C at position 6
    pyrimidine_base_pattern = Chem.MolFromSmarts("n1cncc[c,n]1")
    if not mol.HasSubstructMatch(pyrimidine_base_pattern):
       return False, "No pyrimidine base detected"

    # Specific pyrimidine bases

    cytosine_pattern = Chem.MolFromSmarts("n1ccnc(N)c1=O")
    thymine_pattern = Chem.MolFromSmarts("n1cnc(C)c1=O")
    uracil_pattern = Chem.MolFromSmarts("n1ccnc(O)c1=O")
    
    if not (mol.HasSubstructMatch(cytosine_pattern) or mol.HasSubstructMatch(thymine_pattern) or mol.HasSubstructMatch(uracil_pattern)):
        return False, "No known pyrimidine base found"


    # Deoxyribose sugar, C1 must be linked to N of pyrimidine, C2 must not have an OH group attached
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1([C][C@H]([C@@H](CO)O)O)[C][N]") # 1' is linked to a C (can be anything), no OH at 2', 3'- and 4'-OH, C5'-CH2OH and N for glycosidic bond
    
    if not mol.HasSubstructMatch(deoxyribose_pattern):
         return False, "No deoxyribose sugar detected"


    # Check if there is ribose sugar instead of deoxyribose - fails if this pattern matches,
    ribose_pattern = Chem.MolFromSmarts("[C@H]1([C][C@H]([C@@H](CO)O)[C@@H](O)O)[C][N]") # 1' is linked to a C (can be anything), 2'-OH present, 3'- and 4'-OH, C5'-CH2OH and N for glycosidic bond
    if mol.HasSubstructMatch(ribose_pattern):
        return False, "Ribose detected, not deoxyribose"

    # Check if phosphate is present (indicates nucleotide, not nucleoside), must NOT have P atom
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O])([O])")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group detected, molecule is a nucleotide, not a nucleoside"

    # Verify glycosidic bond (N of base to C1' of sugar)
    glycosidic_bond_pattern = Chem.MolFromSmarts("n~[C@H]1([C][C@H]([C@@H](CO)O)O)[C]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond detected between pyrimidine and sugar"

    return True, "Pyrimidine deoxyribonucleoside detected"
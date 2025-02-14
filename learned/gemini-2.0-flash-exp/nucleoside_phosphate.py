"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:24676 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate consists of a nucleobase, a sugar, and one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for phosphate groups (at least one)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
         return False, "No phosphate group found"

    # 2. Check for a sugar ring (ribose or deoxyribose or other 5-membered ring with 2 or more OH groups)
    sugar_pattern1 = Chem.MolFromSmarts("C1[CH](O)[CH](O)[CH]([CH2X4,CHX4]O)O1") #ribose, deoxyribose or something similar
    sugar_pattern2 = Chem.MolFromSmarts("C1[CH](O)[CH](O)[CH](C)O1") # another potential ribose-like sugar
    sugar_pattern3 = Chem.MolFromSmarts("C1[CH](O)[CH](O)[CH](C(C)=O)O1") # another potential ribose-like sugar

    if not (mol.HasSubstructMatch(sugar_pattern1) or mol.HasSubstructMatch(sugar_pattern2) or mol.HasSubstructMatch(sugar_pattern3)):
       return False, "No sugar ring found"

    # 3. Check for a nucleobase (pyrimidine or purine)
    # Pyrimidines (cytosine, thymine, uracil)
    pyrimidine_pattern1 = Chem.MolFromSmarts("n1cncc(=O)[nH]1") #cytosine/uracil
    pyrimidine_pattern2 = Chem.MolFromSmarts("n1cc(C)c(=O)[nH]1") #thymine
    # Purines (adenine, guanine, hypoxanthine)
    purine_pattern1 = Chem.MolFromSmarts("n1cnc2c1ncnc2N") #adenine
    purine_pattern2 = Chem.MolFromSmarts("n1cnc2c1nc(=O)n[nH]2") #guanine
    purine_pattern3 = Chem.MolFromSmarts("n1cnc2c1ncnc2O") #hypoxanthine

    if not (mol.HasSubstructMatch(pyrimidine_pattern1) or mol.HasSubstructMatch(pyrimidine_pattern2) or
            mol.HasSubstructMatch(purine_pattern1) or mol.HasSubstructMatch(purine_pattern2) or
            mol.HasSubstructMatch(purine_pattern3)):
      return False, "No nucleobase found"

    # 4. Check that sugar is connected to a phosphate and a nucleobase (difficult without complex substructure matching)

    return True, "Contains a nucleobase, a sugar, and at least one phosphate group"
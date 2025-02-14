"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is an iridoid monoterpenoid,
                         (False, reason) otherwise.
                         (None, None) if error.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None #Invalid input
    
    # Basic Iridoid core structure (cyclopentane fused with pyran)
    iridoid_core_pattern1 = Chem.MolFromSmarts("[CH2X4][CHX4]1[CH2X4][OX2][CHX4]2[CHX4]1[CH2X4][CH2X4]2") # basic fused ring
    iridoid_core_pattern2 = Chem.MolFromSmarts("[CHX4]1[CHX4][CHX4][OX2][CHX4]2[CHX4]1[CH2X4][CH2X4]2") #Allow CH at the ring junction
    iridoid_core_pattern3 = Chem.MolFromSmarts("[CH2X4][CHX4]1[CH2X4][OX2][CHX4]2[CHX4]1[CHX4]=[CHX4]2") # one double bond on ring
    iridoid_core_pattern4 = Chem.MolFromSmarts("[CHX4]1[CHX4][CHX4][OX2][CHX4]2[CHX4]1[CHX4]=[CHX4]2") # one double bond on ring
    # Secoiridoid (cleaved cyclopentane ring, pyran) - some flexibility in cleavage location
    secoiridoid_core_pattern1 = Chem.MolFromSmarts("[CH2X4]1[CH2X4][CHX4][OX2][CHX4]2[CHX4]1[CHX4]2") #basic seco, 1 bond broken
    secoiridoid_core_pattern2 = Chem.MolFromSmarts("[CH2X4]1[CHX4][CHX4][OX2][CHX4]2[CHX4]1[CHX4]2") #basic seco, 1 bond broken
    secoiridoid_core_pattern3 = Chem.MolFromSmarts("[CH2X4]1[CHX4][CHX4][OX2][CHX4]2[CHX4]1[CHX4]=[CHX4]2") #basic seco, 1 bond broken, one double bond
    secoiridoid_core_pattern4 = Chem.MolFromSmarts("[CH2X4]1[CHX4]=[CHX4][OX2][CHX4]2[CHX4]1[CHX4]2") #basic seco, 1 bond broken, one double bond

    if not (mol.HasSubstructMatch(iridoid_core_pattern1) or
            mol.HasSubstructMatch(iridoid_core_pattern2) or
            mol.HasSubstructMatch(iridoid_core_pattern3) or
            mol.HasSubstructMatch(iridoid_core_pattern4) or
            mol.HasSubstructMatch(secoiridoid_core_pattern1) or
            mol.HasSubstructMatch(secoiridoid_core_pattern2) or
            mol.HasSubstructMatch(secoiridoid_core_pattern3) or
            mol.HasSubstructMatch(secoiridoid_core_pattern4) ):
        return False, "No core iridoid or secoiridoid structure found"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for carbon count of 10 (ignoring glycosidic parts)
    # Need to try to exclude carbohydrates from the count.
    # The most common glycosidic linkages should be removed from the count.
    # Glucose is a good starting point
    glucose_pattern = Chem.MolFromSmarts("C[C@H]1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    n_glucose_carbons = len(glucose_matches) * 6
    modified_c_count = c_count - n_glucose_carbons

    if modified_c_count < 9 or modified_c_count > 14:  # 10 C +/- a couple of carbons for extra methyls/etc.
         return False, f"Incorrect carbon count {modified_c_count} (glycoside C removed)"
    
    
    return True, "Iridoid or secoiridoid core structure with correct carbon count found"
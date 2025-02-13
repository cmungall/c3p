"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for deoxyribose sugar pattern
    # [OH0,1] means oxygen with 0 or 1 hydrogens (to account for modifications)
    # The pattern looks for a 5-membered ring with an oxygen and appropriate substitutions
    deoxyribose_pattern = Chem.MolFromSmarts("[OH0,1]1[CH2][CH1,2][CH1]([CH2][OH0,1])[CH]1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"

    # Check for pyrimidine base pattern
    # Looking for the basic 6-membered ring with nitrogens at positions 1 and 3
    # This will match cytosine, thymine, uracil and their derivatives
    pyrimidine_pattern = Chem.MolFromSmarts("n1c([*,H])c([*,H])c([*,H])nc1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # Check for N-glycosidic bond between sugar and base
    # The pattern looks for the connection between the sugar and the pyrimidine
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OH0,1]1[CH2][CH1,2][CH1][CH]1[NH0]2[c]([*,H])[c]([*,H])[c]([*,H])[n][c]2")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No N-glycosidic bond between sugar and pyrimidine found"

    # Additional check for carbonyl groups in the pyrimidine ring
    # Pyrimidine bases should have at least one C=O group
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing characteristic carbonyl group in pyrimidine base"

    # Check ring count to avoid complex molecules with additional large rings
    rings = mol.GetRingInfo().NumRings()
    if rings > 3:  # Should typically have 2 rings (sugar + base)
        # Allow up to 3 rings to account for some modifications
        if rings > 4:
            return False, "Too many rings in structure"

    return True, "Contains deoxyribose sugar connected to pyrimidine base via N-glycosidic bond"
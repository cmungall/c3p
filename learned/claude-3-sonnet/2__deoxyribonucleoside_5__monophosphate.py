"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:73333 2'-deoxyribonucleoside 5'-monophosphate
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
    
    # Look for deoxyribose pattern (furanose ring with missing 2' OH)
    # [CH2]-[CH]-[CH](-[OH])-[CH]-O (cyclic)
    deoxyribose_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH]([OX2H1])[CH]1O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"
    
    # Look for phosphate group (-O-P(=O)(-O)(-O)-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatches(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for nucleobase attachment via N-glycosidic bond
    # Pattern matches N connected to sugar carbon
    n_glycosidic_pattern = Chem.MolFromSmarts("[NX3]1@[CH]2@[CH]@[CH]@[CH]@O2")
    if not mol.HasSubstructMatch(n_glycosidic_pattern):
        return False, "No N-glycosidic bond found"
    
    # Verify phosphate is at 5' position by checking connection to CH2
    phosphate_5_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[PX4](=[OX1])([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(phosphate_5_pattern):
        return False, "Phosphate not at 5' position"
    
    # Look for common nucleobase patterns (purine or pyrimidine)
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]cnc12")  # adenine-like
    pyrimidine_pattern = Chem.MolFromSmarts("c1cn([*])c(=O)[nH]c1=O")  # thymine/uracil-like
    cytosine_pattern = Chem.MolFromSmarts("Nc1cc[nH]c(=O)n1")  # cytosine-like
    
    if not (mol.HasSubstructMatch(purine_pattern) or 
            mol.HasSubstructMatch(pyrimidine_pattern) or 
            mol.HasSubstructMatch(cytosine_pattern)):
        return False, "No recognizable nucleobase found"
    
    # Count phosphorus atoms - should only have one
    p_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15])
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, should be exactly 1"
    
    return True, "Contains deoxyribose sugar with 5' phosphate and nucleobase connected via N-glycosidic bond"
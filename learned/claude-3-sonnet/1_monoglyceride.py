"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at position 1
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has:
    - A glycerol backbone (3 carbons)
    - An ester group at position 1
    - Hydroxyl groups at positions 2 and 3
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_1_monoglyceride, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for 1-monoglyceride:
    # [#6]-C(=O)-O-CH2-CH(OH)-CH2OH
    # This pattern captures the essential features:
    # - Ester group at position 1
    # - Hydroxyl groups at positions 2 and 3
    # - Carbon chain attached to the carbonyl
    pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[CH2X4]-[CH1X4](-[OX2H1])-[CH2X4]-[OX2H1]")
    
    if pattern is None:
        return False, "Invalid SMARTS pattern"
        
    if not mol.HasSubstructMatch(pattern):
        return False, "Does not match 1-monoglyceride pattern"

    # Count ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if ester_pattern is None:
        return False, "Invalid ester SMARTS pattern"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    if hydroxyl_pattern is None:
        return False, "Invalid hydroxyl SMARTS pattern"
    
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need exactly 2"

    # Verify carbon chain length (should be at least 2 carbons in acyl group)
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[#6]-[#6]")
    if acyl_pattern is None:
        return False, "Invalid acyl SMARTS pattern"
        
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with one ester bond at position 1 and two free hydroxyl groups"
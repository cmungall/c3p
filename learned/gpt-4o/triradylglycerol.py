"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Pattern to identify ester bond
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    # Patterns for acyl, alkyl, and alk-1-enyl groups
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    alkyl_pattern = Chem.MolFromSmarts("[C][O]")  # Extended to match different possible alkyl groups
    alk1enyl_pattern = Chem.MolFromSmarts("[C]=[C][O]")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    alk1enyl_matches = mol.GetSubstructMatches(alk1enyl_pattern)

    # Check for at least three ester or possible linkage groups
    if len(ester_matches) < 3:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 3"

    # Ensure diverse presence of substituents across positions
    if len(acyl_matches) > 0 or len(alkyl_matches) > 0 or len(alk1enyl_matches) > 0:
        if len(acyl_matches) + len(alkyl_matches) + len(alk1enyl_matches) >= 3:
            return True, "Is a triradylglycerol with enough diversity in substituents"
    
    return False, "Insufficient diversity among acyl/alkyl/alk-1-enyl substituents"
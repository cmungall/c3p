"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:57980 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create and verify SMARTS patterns
    patterns = {
        'phosphate': '[OX2,O-]P(=O)([OX2,O-])[OX2]',  # Phosphate group
        'ester': '[CX3](=O)[OX2]',  # Ester group
        'glycerol_backbone': '[CH2X4,CH2][CH,CHX4][CH2X4,CH2]',  # More flexible glycerol pattern
        'ethanolamine': 'OCCN',  # Pattern to exclude
        'choline': 'OCC[N+]',    # Pattern to exclude
        'ether': '[CH2,CH][OX2][CH2]C'  # Pattern to exclude ether lipids
    }
    
    compiled_patterns = {}
    for name, pattern in patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat is None:
            return False, f"Invalid SMARTS pattern for {name}"
        compiled_patterns[name] = pat

    # Exclude molecules with ethanolamine or choline groups
    if mol.HasSubstructMatch(compiled_patterns['ethanolamine']):
        return False, "Contains ethanolamine group"
    if mol.HasSubstructMatch(compiled_patterns['choline']):
        return False, "Contains choline group"

    # Check for phosphate group
    if not mol.HasSubstructMatch(compiled_patterns['phosphate']):
        return False, "No phosphate group found"

    # Check for single ester group (acyl chain)
    ester_matches = mol.GetSubstructMatches(compiled_patterns['ester'])
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(compiled_patterns['glycerol_backbone']):
        return False, "No glycerol backbone found"

    # Exclude ether lipids
    if mol.HasSubstructMatch(compiled_patterns['ether']):
        return False, "Contains ether linkage"

    # Count oxygen atoms - should be 7 
    # (3 from phosphate, 2 from ester, 2 from glycerol backbone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 7:
        return False, f"Found {o_count} oxygen atoms, expected 7"

    # Count phosphorus atoms - should be exactly 1
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, need exactly 1"

    # Additional check for acyl chain length (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:  # minimum 3 for glycerol + 4 for shortest acyl chain
        return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with one acyl chain and phosphate group"
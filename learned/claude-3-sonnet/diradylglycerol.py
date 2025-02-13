"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:18035 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups (acyl, alkyl, or alk-1-enyl)
    at any two of the three possible positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Look for 2 substituent groups (acyl, alkyl, or alk-1-enyl)
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    alkyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkenyl_pattern = Chem.MolFromSmarts("[CX3]=[CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    
    substituent_matches = mol.GetSubstructMatches(acyl_pattern) + mol.GetSubstructMatches(alkyl_pattern) + mol.GetSubstructMatches(alkenyl_pattern)
    if len(substituent_matches) < 2:
        return False, f"Found {len(substituent_matches)} substituent groups, need at least 2"

    # Check molecular weight - diradylglycerols typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for diradylglycerol"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for diradylglycerol"
    if o_count != 5:
        return False, "Must have exactly 5 oxygens (2 ester groups, 1 free hydroxyl)"

    return True, "Contains glycerol backbone with 2 substituent groups at any 2 positions"
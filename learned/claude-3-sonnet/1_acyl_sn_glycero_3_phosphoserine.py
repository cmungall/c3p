"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:75944 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with specific substitution pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphoserine group
    # [P](=O)(O)-O-C-C(N)-C(=O)O
    phosphoserine_pattern = Chem.MolFromSmarts("[PX4](=O)([OH])[OX2][CH2][CH]([NH2])[CX3](=O)[OH]")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Check for single ester group (acyl chain)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for free hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[CHX4][OHX2]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No free hydroxyl group found"

    # Check for fatty acid chain (at least 4 carbons in chain)
    fatty_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "No fatty acid chain found"

    # Count key atoms to verify overall composition
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if n_count != 1:
        return False, "Must have exactly one nitrogen atom"

    # Verify stereochemistry if specified
    if '@' in smiles:
        # The glycerol carbon should be R configuration (@@ in SMILES) if specified
        # and the serine carbon should be S configuration (@ in SMILES) if specified
        if not ('@@' in smiles and '@' in smiles and '@@' in smiles.split('@')[1]):
            return False, "Incorrect stereochemistry for 1-acyl-sn-glycero-3-phosphoserine"

    return True, "Contains glycerol backbone with single acyl chain and phosphoserine group in correct configuration"
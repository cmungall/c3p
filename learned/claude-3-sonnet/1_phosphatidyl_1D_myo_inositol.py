"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:16813 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for myo-inositol ring
    # Pattern matches cyclohexane ring with 6 OH groups in myo configuration
    myo_inositol_pattern = Chem.MolFromSmarts("[C@@H]1([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]1[OH]")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol ring found or incorrect stereochemistry"

    # Check for phosphate group (-P(=O)(O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([O,OH])[O]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone with correct stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found or incorrect stereochemistry"

    # Check for two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Verify overall composition
    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if o_count < 12:  # 6 from inositol, 2 from esters, 4 from phosphate
        return False, "Insufficient oxygen atoms"
    if c_count < 25:  # minimum for shortest possible chains
        return False, "Carbon chains too short"

    # Check molecular weight - these are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low"

    return True, "Contains myo-inositol ring connected via phosphate to glycerol with two fatty acid chains"
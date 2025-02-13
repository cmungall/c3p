"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: CHEBI:18294 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerolipid with a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(~O)(~O)(~O)(~O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for ester linkage between phosphate and glycerol
    ester_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]P(~O)(~O)(~O)(~O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Phosphate not ester-linked to glycerol backbone"
    
    # Look for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"
    
    # Check molecular weight - glycerophospholipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycerophospholipid"
    
    # Count carbons, oxygens and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20:
        return False, "Too few carbons for glycerophospholipid"
    if o_count < 6:
        return False, "Too few oxygens for glycerophospholipid"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus atom"

    return True, "Contains glycerol backbone with phosphate group and fatty acid chains"
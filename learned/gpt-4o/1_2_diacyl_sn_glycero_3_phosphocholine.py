"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (central CH with two O and rest atoms connected)
    glycerol_backbone_pattern = Chem.MolFromSmarts("[C@H](COC(=O)*)OC(=O)*")
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphocholine group: [O-]P(=O)(OCC[N+](C)(C)C)
    phosphocholine_group_pattern = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)")
    if not mol.HasSubstructMatch(phosphocholine_group_pattern):
        return False, "No phosphocholine group found"
    
    # Check for negative charge balance
    net_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if net_charge != 0:
        return False, f"Net charge imbalance: {net_charge} (expected 0)"
    
    # Check the presence of exactly two ester linkages
    ester_linkage_pattern = Chem.MolFromSmarts("COC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    return True, "Contains 1,2-diacyl-sn-glycero-3-phosphocholine features"
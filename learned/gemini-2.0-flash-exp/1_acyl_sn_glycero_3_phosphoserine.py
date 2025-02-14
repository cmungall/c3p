"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoserine has a glycerol backbone with an acyl group at the 1-position,
    and a phosphoserine group at the 3-position.

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

    # Check for glycerol backbone with sn-1 and sn-3 stereochemistry
    # Specifically looking for stereochem. of -C[C@H](O)C- 
    glycerol_pattern = Chem.MolFromSmarts("[CX4][C@H](O)[CX4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-glycerol backbone found"
    
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    glycerol_c2_idx = glycerol_match[1] # Get the index of the center carbon of glycerol
    glycerol_c1_idx = glycerol_match[0] # Get the index of the first carbon of glycerol
    glycerol_c3_idx = glycerol_match[2] # Get the index of the third carbon of glycerol
    
    # Check if carbons are connected to each other and in order
    if not mol.GetBondBetweenAtoms(glycerol_c1_idx,glycerol_c2_idx):
        return False, "Glycerol atoms are not connected correctly"
    if not mol.GetBondBetweenAtoms(glycerol_c2_idx,glycerol_c3_idx):
        return False, "Glycerol atoms are not connected correctly"


    # Check for ester group attached to sn-1 carbon
    ester_pattern = Chem.MolFromSmarts("[CX4](-[OX2])=[OX1]")
    
    
    sn1_ester_match = None
    for bond in mol.GetAtomWithIdx(glycerol_c1_idx).GetBonds():
        if bond.GetOtherAtomIdx(glycerol_c1_idx) != glycerol_c2_idx:
            if mol.GetAtomWithIdx(bond.GetOtherAtomIdx(glycerol_c1_idx)).HasSubstructMatch(ester_pattern):
                sn1_ester_match = mol.GetAtomWithIdx(bond.GetOtherAtomIdx(glycerol_c1_idx)).GetSubstructMatch(ester_pattern)
                break
    if sn1_ester_match is None:
       return False, "No ester group at sn-1 position"

    
    # check if chain connected to ester is long enough (>=4 carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    
    acyl_chain_matches = None
    for atom in mol.GetAtomWithIdx(sn1_ester_match[0]).GetNeighbors():
        if atom.GetIdx() != glycerol_c1_idx:
            if mol.GetSubstructMatch(acyl_chain_pattern, atom.GetIdx()) is not None:
              acyl_chain_matches = mol.GetSubstructMatch(acyl_chain_pattern,atom.GetIdx())
              break;
    if acyl_chain_matches is None:
        return False, "No acyl chain connected to the ester at sn-1 position or chain is too short"



    # Check for phosphoserine group at sn-3 carbon
    phosphoserine_pattern = Chem.MolFromSmarts("P(=O)(O)-O[CX4][CX4](N)C(=O)O")
    phosphoserine_match = None;
    for bond in mol.GetAtomWithIdx(glycerol_c3_idx).GetBonds():
         if bond.GetOtherAtomIdx(glycerol_c3_idx) != glycerol_c2_idx:
            if mol.GetAtomWithIdx(bond.GetOtherAtomIdx(glycerol_c3_idx)).HasSubstructMatch(phosphoserine_pattern):
               phosphoserine_match = mol.GetAtomWithIdx(bond.GetOtherAtomIdx(glycerol_c3_idx)).GetSubstructMatch(phosphoserine_pattern)
               break
    if phosphoserine_match is None:
      return False, "No phosphoserine group found at sn-3 position"



    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if n_count != 1:
         return False, f"Must have exactly 1 nitrogen, got {n_count}"

    if p_count != 1:
         return False, f"Must have exactly 1 phosphorus, got {p_count}"
    if o_count < 7:
        return False, f"Must have at least 7 oxygens for 3-phosphoserine, got {o_count}"

    if c_count < 6:
        return False, f"Too few carbons, got {c_count}"


    return True, "Matches 1-acyl-sn-glycero-3-phosphoserine criteria"
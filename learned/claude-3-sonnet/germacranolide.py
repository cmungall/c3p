"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: germacranolide
A sesquiterpene lactone based on germacrane skeleton
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for lactone ring (furanone)
    lactone_pattern = Chem.MolFromSmarts("[C]1[C](=[O])[O][C]1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for 10-membered ring (germacrane skeleton)
    # Look for connected atoms forming a 10-membered ring
    ring_info = mol.GetRingInfo()
    has_10_ring = False
    for ring in ring_info.AtomRings():
        if len(ring) == 10:
            has_10_ring = True
            break
    if not has_10_ring:
        return False, "No 10-membered ring found"

    # Count carbons (should be around 15 for sesquiterpene)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 17:  # Allow some flexibility for modifications
        return False, f"Carbon count {c_count} not consistent with sesquiterpene lactone"

    # Check for α-methylene-γ-lactone pattern (common in germacranolides)
    # but not required as some germacranolides might have variations
    methylene_lactone = Chem.MolFromSmarts("[CH2]=[C]1[C](=[O])[O][C]1")
    has_methylene = mol.HasSubstructMatch(methylene_lactone)

    # Check for ring fusion between lactone and 10-membered ring
    fused_ring_pattern = Chem.MolFromSmarts("[C]1[C](=[O])[O][C]1[C]2[C,c][C,c][C,c][C,c][C,c][C,c][C,c][C,c][C,c]2")
    has_fusion = mol.HasSubstructMatch(fused_ring_pattern)
    if not has_fusion:
        return False, "Lactone not properly fused to 10-membered ring"

    # Calculate number of rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 2:  # Should have at least 2 rings (10-membered + lactone)
        return False, "Insufficient number of rings"

    # Additional check for oxygen count (should have at least 2 for lactone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for germacranolide"

    # If we've made it here, we have the basic germacranolide structure
    reason = "Contains germacrane skeleton (10-membered ring) fused to lactone"
    if has_methylene:
        reason += " with α-methylene-γ-lactone"
    
    return True, reason
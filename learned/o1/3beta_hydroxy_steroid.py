"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36804 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group at the 3-position in the beta orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone: 4 fused rings (three 6-membered and one 5-membered)
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    if len(rings) < 4:
        return False, "Less than 4 rings found"
    
    # Count ring sizes
    ring_sizes = [len(r) for r in rings]
    count_6_membered = ring_sizes.count(6)
    count_5_membered = ring_sizes.count(5)
    if count_6_membered < 3 or count_5_membered < 1:
        return False, "Does not have three 6-membered rings and one 5-membered ring"
    
    # Check that the rings are fused together
    fused = True
    for i, ring1 in enumerate(rings):
        for ring2 in rings[i+1:]:
            shared_atoms = set(ring1) & set(ring2)
            if len(shared_atoms) == 0:
                fused = False
                break
        if not fused:
            break
    if not fused:
        return False, "The rings are not fully fused together"
    
    # Look for 3beta-hydroxy group
    # In steroids, C3 is typically in ring A (one of the six-membered rings)
    # We'll look for a chiral carbon connected to an OH group in the ring system
    pattern = Chem.MolFromSmarts('[C@H](O)[#6]')
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No chiral carbon with OH group found"
    
    # Verify that the chiral carbon is part of the ring system
    found_beta_oh = False
    for match in matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        # Check if the atom is in the ring system
        if any(atom_idx in ring for ring in rings):
            # Check chiral tag for beta orientation
            chiral_tag = atom.GetChiralTag()
            if chiral_tag == Chem.CHI_TETRAHEDRAL_CCW:
                found_beta_oh = True
                break
    if not found_beta_oh:
        return False, "No 3beta-hydroxy group found"
    
    return True, "Contains steroid backbone with 3beta-hydroxy group"
"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify sugar rings (5 or 6-membered rings with one oxygen atom)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_rings = []
    for ring in rings:
        if len(ring) not in [5, 6]:
            continue  # Not a 5 or 6-membered ring
        o_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                o_count += 1
            elif atomic_num == 6:
                c_count += 1
            else:
                break  # Contains other atoms, not a sugar ring
        else:
            if o_count == 1 and c_count == len(ring) - 1:
                sugar_rings.append(set(ring))

    if len(sugar_rings) < 2:
        return False, f"Found {len(sugar_rings)} sugar rings, need at least 2"

    # Map atom indices to sugar ring indices
    ring_atom_map = {}
    for ring_idx, ring_atoms in enumerate(sugar_rings):
        for atom_idx in ring_atoms:
            ring_atom_map[atom_idx] = ring_idx

    # Find glycosidic bonds connecting two different sugar rings
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        idx1 = atom1.GetIdx()
        idx2 = atom2.GetIdx()
        atomic_num1 = atom1.GetAtomicNum()
        atomic_num2 = atom2.GetAtomicNum()

        # Check for C-O bonds
        if {atomic_num1, atomic_num2} != {6, 8}:
            continue  # Not a C-O bond

        # Identify carbon and oxygen atoms
        if atomic_num1 == 8:
            oxygen_atom = atom1
            carbon_atom = atom2
        else:
            oxygen_atom = atom2
            carbon_atom = atom1

        # Oxygen atom should not be in a ring (glycosidic oxygen)
        if ring_info.IsAtomInRingOfSize(oxygen_atom.GetIdx(), 5) or ring_info.IsAtomInRingOfSize(oxygen_atom.GetIdx(), 6):
            continue

        # Carbon atom should be in a ring (part of a sugar ring)
        if carbon_atom.GetIdx() not in ring_atom_map:
            continue

        # Find the other atom connected to the oxygen
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() == carbon_atom.GetIdx():
                continue  # Skip the already considered carbon atom
            if neighbor.GetAtomicNum() != 6:
                continue  # Not connected to another carbon atom

            # Neighboring carbon atom should be in a ring
            if neighbor.GetIdx() not in ring_atom_map:
                continue

            # Ensure the two carbons are from different sugar rings
            ring_idx1 = ring_atom_map[carbon_atom.GetIdx()]
            ring_idx2 = ring_atom_map[neighbor.GetIdx()]
            if ring_idx1 != ring_idx2:
                glycosidic_bonds.append((oxygen_atom.GetIdx(), carbon_atom.GetIdx(), neighbor.GetIdx()))

    if len(glycosidic_bonds) == 0:
        return False, "No glycosidic bond connecting sugar rings found"

    # Ensure only two sugar rings are present
    connected_rings = set()
    for _, c_idx1, c_idx2 in glycosidic_bonds:
        ring_idx1 = ring_atom_map[c_idx1]
        ring_idx2 = ring_atom_map[c_idx2]
        connected_rings.update([ring_idx1, ring_idx2])

    if len(connected_rings) != 2 or len(sugar_rings) > 2:
        return False, "More than two sugar rings connected; not a disaccharide"

    return True, "Molecule is a disaccharide with two monosaccharide units connected via a glycosidic bond"
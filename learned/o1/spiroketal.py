"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Build a dictionary to count the number of rings each atom is in
    atom_ring_count = {}
    for ring in atom_rings:
        for atom_idx in ring:
            atom_ring_count[atom_idx] = atom_ring_count.get(atom_idx, 0) + 1

    # Identify potential spiro centers (atoms in exactly two rings)
    spiro_atoms = [atom_idx for atom_idx, count in atom_ring_count.items() if count == 2]

    if not spiro_atoms:
        return False, "No spiro centers found"

    # Loop over spiro atoms to check for ketal functionality
    for atom_idx in spiro_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)

        # Check if the two rings share only this atom
        rings_with_atom = [ring for ring in atom_rings if atom_idx in ring]
        if len(rings_with_atom) != 2:
            continue

        shared_atoms = set(rings_with_atom[0]) & set(rings_with_atom[1])
        if len(shared_atoms) != 1:
            continue  # The rings share more than one atom

        # Check if the atom is carbon
        if atom.GetAtomicNum() != 6:
            continue  # Not a carbon atom

        # Check if the atom is connected to exactly two oxygens
        neighbors = atom.GetNeighbors()
        oxygens = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            continue  # Not connected to exactly two oxygens

        # Check that each oxygen is connected to the spiro carbon and another carbon (ether linkage)
        is_ketal = True
        for oxygen in oxygens:
            if oxygen.GetDegree() != 2:
                is_ketal = False
                break
            # Get the other atom connected to oxygen (excluding spiro carbon)
            oxygen_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetIdx() != atom_idx]
            if len(oxygen_neighbors) != 1:
                is_ketal = False
                break
            if oxygen_neighbors[0].GetAtomicNum() != 6:
                is_ketal = False
                break
        if not is_ketal:
            continue  # Does not fulfill ketal criteria

        # If all checks passed, this is a spiroketal
        return True, f"Molecule is a spiroketal with spiro ketal center at atom index {atom_idx}"

    return False, "No spiroketal functionality found"
"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
"""

from rdkit import Chem

def get_fused_ring_subgraphs(mol):
    """
    Helper function to identify fused ring systems in a molecule.

    Args:
        mol: rdkit.Chem.Mol object

    Returns:
        List of sets, each set contains indices of rings that are fused together
    """
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Build adjacency lists for rings
    ring_adj_list = []
    num_rings = len(atom_rings)
    for i in range(num_rings):
        ring_i_atoms = set(atom_rings[i])
        fused = []
        for j in range(num_rings):
            if i != j:
                ring_j_atoms = set(atom_rings[j])
                if ring_i_atoms & ring_j_atoms:
                    fused.append(j)
        ring_adj_list.append(fused)

    # Find connected components in ring adjacency graph
    from collections import deque
    visited = set()
    fused_ring_systems = []
    for i in range(num_rings):
        if i not in visited:
            system = set()
            queue = deque([i])
            while queue:
                idx = queue.popleft()
                if idx not in visited:
                    visited.add(idx)
                    system.add(idx)
                    neighbors = ring_adj_list[idx]
                    for n in neighbors:
                        if n not in visited:
                            queue.append(n)
            fused_ring_systems.append(system)
    return fused_ring_systems

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a steroid molecule where at least one hydroxy group is esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Get fused ring systems
    fused_ring_systems = get_fused_ring_subgraphs(mol)

    # Identify steroid backbone
    steroid_found = False
    steroid_atom_indices = set()

    for ring_idxs in fused_ring_systems:
        # Collect atom indices and ring sizes
        ring_atoms = set()
        ring_sizes = []
        for ring_idx in ring_idxs:
            ring = atom_rings[ring_idx]
            ring_atoms.update(ring)
            ring_sizes.append(len(ring))
        # Check for 4-ring system with sizes 6,6,6,5 and 17 carbons
        if len(ring_idxs) == 4:
            ring_sizes_sorted = sorted(ring_sizes)
            if ring_sizes_sorted == [5,6,6,6] and len(ring_atoms) == 17:
                steroid_found = True
                steroid_atom_indices = ring_atoms
                break

    if not steroid_found:
        return False, "No steroid backbone found"

    # Define sulfate ester group patterns (protonated and deprotonated)
    sulfate_pattern = Chem.MolFromSmarts('OS(=O)(=O)[O-]')
    sulfate_pattern_unprotonated = Chem.MolFromSmarts('OS(=O)(=O)O')

    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern) + mol.GetSubstructMatches(sulfate_pattern_unprotonated)

    if not sulfate_matches:
        return False, "No sulfate ester group found"

    # Check if sulfate is attached to steroid backbone via oxygen
    for match in sulfate_matches:
        ester_oxygen_idx = match[0]  # Oxygen connected to sulfur
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        # Check neighbors of ester oxygen atom
        for neighbor in ester_oxygen_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in steroid_atom_indices:
                return True, "Contains steroid backbone with sulfate ester group attached via oxygen"

    return False, "Sulfate group not attached to steroid backbone via ester linkage"
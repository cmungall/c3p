"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is based on the cyclopenta[a]phenanthrene skeleton,
    consisting of three fused six-membered rings and one fused five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        return False, "No rings found in the molecule"

    # Build a ring adjacency graph
    import collections
    ring_graph = collections.defaultdict(set)
    for i, ring1 in enumerate(atom_rings):
        for j, ring2 in enumerate(atom_rings):
            if i >= j:
                continue
            # Check if rings share at least two atoms (fused rings)
            shared_atoms = set(ring1) & set(ring2)
            if len(shared_atoms) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph (fused ring systems)
    def dfs(node, visited, component):
        visited.add(node)
        component.add(node)
        for neighbor in ring_graph[node]:
            if neighbor not in visited:
                dfs(neighbor, visited, component)
    
    visited = set()
    fused_ring_systems = []
    for node in ring_graph:
        if node not in visited:
            component = set()
            dfs(node, visited, component)
            fused_ring_systems.append(component)
    
    # Look for a fused ring system with four rings of sizes 6,6,6,5
    found_steroid_core = False
    for system in fused_ring_systems:
        if len(system) != 4:
            continue
        ring_sizes = [len(atom_rings[i]) for i in system]
        # Check for three six-membered rings and one five-membered ring
        if sorted(ring_sizes) == [5,6,6,6]:
            found_steroid_core = True
            break
    
    if not found_steroid_core:
        return False, "Steroid core not found (requires fused rings of sizes 6,6,6,5)"

    return True, "Contains steroid core (cyclopenta[a]phenanthrene skeleton)"
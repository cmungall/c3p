"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for triterpenoid core
    # Triterpenoids usually have 4 or more fused rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Contains {num_rings} rings, less than 4 rings required for triterpenoid core"

    # Identify fused ring systems
    rings = ring_info.AtomRings()
    ring_sets = [set(ring) for ring in rings]
    ring_adjacency = {i: set() for i in range(len(ring_sets))}
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            if len(ring_sets[i].intersection(ring_sets[j])) >= 2:
                ring_adjacency[i].add(j)
                ring_adjacency[j].add(i)

    def get_fused_ring_systems(ring_adjacency):
        visited = set()
        fused_ring_systems = []
        for ring_idx in ring_adjacency:
            if ring_idx not in visited:
                # Start a new fused ring system
                stack = [ring_idx]
                fused_system = set()
                while stack:
                    curr = stack.pop()
                    if curr not in visited:
                        visited.add(curr)
                        fused_system.add(curr)
                        stack.extend(ring_adjacency[curr] - visited)
                fused_ring_systems.append(fused_system)
        return fused_ring_systems

    fused_ring_systems = get_fused_ring_systems(ring_adjacency)

    # Get the largest fused ring system
    largest_fused_system = max(fused_ring_systems, key=lambda x: len(x))
    num_fused_rings = len(largest_fused_system)

    if num_fused_rings < 4:
        return False, f"Insufficient fused rings for triterpenoid core, found only {num_fused_rings} fused rings"

    # Check for sugar moieties (rings of size 5 or 6 containing oxygen)
    sugar_rings_found = False
    for ring in rings:
        ring_size = len(ring)
        if ring_size == 5 or ring_size == 6:
            num_oxygen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() ==8)
            if num_oxygen >= 1:
                sugar_rings_found = True
                break

    if not sugar_rings_found:
        return False, "No sugar moieties found"

    # Simplified glycosidic bond check
    # Assume that if both triterpenoid core and sugar moiety are present, glycosidic bond exists

    return True, "Contains triterpenoid core with sugar moiety"
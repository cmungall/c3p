"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:26870 triterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene, typically containing 30 carbons,
    but may have modifications such as rearrangements, glycosylations, or additions/removals of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27 or c_count > 50:
        return False, f"Carbon count {c_count} not in typical triterpenoid range (27-50 carbons)"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, likely not a triterpenoid"

    # Calculate number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Only {num_rings} rings found, triterpenoids typically have at least 4 rings"

    # Check for fused ring systems
    # Get ring bonds
    ssr = Chem.GetSymmSSSR(mol)
    fused_rings = [set(ring) for ring in ssr]

    # Build a graph of ring adjacencies
    ring_adjacency = []
    for i in range(len(fused_rings)):
        neighbors = []
        for j in range(len(fused_rings)):
            if i != j and len(fused_rings[i] & fused_rings[j]) >= 2:
                neighbors.append(j)
        ring_adjacency.append(neighbors)

    # Find the largest set of fused rings using DFS
    max_fused_rings = 0
    visited = set()
    def dfs(ring_idx):
        visited.add(ring_idx)
        total = 1
        for nbr in ring_adjacency[ring_idx]:
            if nbr not in visited:
                total += dfs(nbr)
        return total

    for i in range(len(fused_rings)):
        if i not in visited:
            fused_ring_count = dfs(i)
            if fused_ring_count > max_fused_rings:
                max_fused_rings = fused_ring_count
                
    if max_fused_rings < 4:
        return False, f"Only {max_fused_rings} fused rings found, triterpenoids typically have at least 4 fused rings"

    # Calculate fraction of sp3 carbons
    sp3_carbons = 0
    total_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            total_carbons += 1
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                sp3_carbons += 1
    sp3_fraction = sp3_carbons / total_carbons if total_carbons > 0 else 0
    if sp3_fraction < 0.6:
        return False, f"Fraction of sp3 carbons is {sp3_fraction:.2f}, which is too low for typical triterpenoids"

    return True, "Molecule meets criteria for a triterpenoid (carbon count, oxygen atoms, fused ring structures, saturation)"
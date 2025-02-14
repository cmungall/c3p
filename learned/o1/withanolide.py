"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as any steroid lactone that is a C28 steroid 
    with a modified side chain forming a lactone ring and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for tetracyclic steroid nucleus (fused ring system with 4 rings: 6-6-6-5)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()
    
    # Build a list of rings and their sizes
    ring_sizes = [len(ring) for ring in atom_rings]
    num_rings = len(atom_rings)
    
    # Build a connectivity map of rings
    ring_adj_matrix = [[False]*num_rings for _ in range(num_rings)]
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # If rings share bonds, they are fused
            if set(bond_rings[i]) & set(bond_rings[j]):
                ring_adj_matrix[i][j] = True
                ring_adj_matrix[j][i] = True

    # Function to find connected components (fused ring systems)
    def get_fused_systems():
        visited = set()
        fused_systems = []
        for i in range(num_rings):
            if i not in visited:
                stack = [i]
                fused_system = []
                while stack:
                    idx = stack.pop()
                    if idx not in visited:
                        visited.add(idx)
                        fused_system.append(idx)
                        for j in range(num_rings):
                            if ring_adj_matrix[idx][j] and j not in visited:
                                stack.append(j)
                fused_systems.append(fused_system)
        return fused_systems

    fused_systems = get_fused_systems()
    has_steroid_nucleus = False
    steroid_atoms = set()

    for system in fused_systems:
        if len(system) == 4:
            system_ring_sizes = [ring_sizes[idx] for idx in system]
            system_ring_sizes.sort()
            if system_ring_sizes == [5,6,6,6]:
                has_steroid_nucleus = True
                # Collect atoms in the steroid nucleus
                for idx in system:
                    steroid_atoms.update(atom_rings[idx])
                break

    if not has_steroid_nucleus:
        return False, "No tetracyclic steroid nucleus (rings of sizes 6-6-6-5) found"

    # Check for lactone ring in side chain (cyclic ester not part of steroid nucleus)
    # Define lactone pattern (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts('C(=O)O[C;R]')
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    has_lactone = False

    for match in lactone_matches:
        lactone_o = match[2]
        # Check if the oxygen is part of a ring
        for ring in atom_rings:
            if lactone_o in ring:
                # Check if ring is not part of steroid nucleus
                ring_atoms_set = set(ring)
                if not ring_atoms_set.issubset(steroid_atoms):
                    has_lactone = True
                    break
        if has_lactone:
            break

    if not has_lactone:
        return False, "No lactone ring found in side chain"

    # Check number of carbons (should be around 28)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 26 or c_count > 30:
        return False, f"Carbon count ({c_count}) not consistent with C28 steroid"

    return True, "Molecule is a withanolide (tetracyclic steroid nucleus with side-chain lactone ring)"
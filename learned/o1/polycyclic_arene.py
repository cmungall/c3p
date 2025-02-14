"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) based on its SMILES string.
    A polycyclic arene is a hydrocarbon consisting of fused aromatic rings, containing only carbon and hydrogen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule contains only carbon and hydrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1:
            return False, "Molecule contains atoms other than carbon and hydrogen"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    num_rings = len(atom_rings)
    if num_rings < 2:
        return False, "Molecule does not have at least two rings"

    # Build ring adjacency graph
    ring_graph = {}
    for i in range(num_rings):
        ring_graph[i] = set()

    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # Check if rings i and j share at least one bond
            # Get bonds in ring i
            bonds_i = set()
            atoms_i = atom_rings[i]
            for idx in range(len(atoms_i)):
                a1 = atoms_i[idx]
                a2 = atoms_i[(idx + 1) % len(atoms_i)]
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond:
                    bonds_i.add(bond.GetIdx())
            # Get bonds in ring j
            bonds_j = set()
            atoms_j = atom_rings[j]
            for idx in range(len(atoms_j)):
                a1 = atoms_j[idx]
                a2 = atoms_j[(idx + 1) % len(atoms_j)]
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond:
                    bonds_j.add(bond.GetIdx())
            # Check for shared bonds
            if bonds_i & bonds_j:
                # Rings i and j share at least one bond (fused)
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # Find connected components (fused ring systems)
    visited = [False] * num_rings
    ring_systems = []
    for i in range(num_rings):
        if not visited[i]:
            stack = [i]
            component = []
            while stack:
                node = stack.pop()
                if not visited[node]:
                    visited[node] = True
                    component.append(node)
                    stack.extend(ring_graph[node])
            ring_systems.append(component)

    # For each ring system, check if it contains at least two aromatic rings made up of carbon atoms
    for system in ring_systems:
        if len(system) < 2:
            continue  # Skip ring systems with less than two rings
        all_rings_are_aromatic = True
        for ring_idx in system:
            atoms_in_ring = atom_rings[ring_idx]
            ring_aromatic = True
            for atom_idx in atoms_in_ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if not atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                    ring_aromatic = False
                    break
            if not ring_aromatic:
                all_rings_are_aromatic = False
                break
        if all_rings_are_aromatic:
            # Found a fused aromatic ring system with at least two rings made up of carbon atoms
            return True, "Molecule is a polycyclic aromatic hydrocarbon"
    
    return False, "No fused aromatic hydrocarbon ring system with at least two aromatic rings found"
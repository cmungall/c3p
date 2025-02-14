"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:28053 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    monosaccharide_rings = []
    for ring in atom_rings:
        if len(ring) == 5 or len(ring) == 6:
            # Check ring composition
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_oxygens_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            num_carbons_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            if num_oxygens_in_ring == 1 and num_carbons_in_ring == len(ring) - 1:
                # Potential monosaccharide ring
                # Check for exocyclic oxygen substituents
                num_exocyclic_oxygens = 0
                exocyclic_oxygens = []
                for atom in ring_atoms:
                    if atom.GetAtomicNum() == 6:
                        # For carbon atoms in the ring
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                                # Exocyclic oxygen found
                                num_exocyclic_oxygens += 1
                                exocyclic_oxygens.append((atom.GetIdx(), neighbor.GetIdx()))
                if num_exocyclic_oxygens >= 3:
                    # Accept this ring as monosaccharide
                    monosaccharide_rings.append({
                        'ring_atoms': set(ring),
                        'exocyclic_oxygens': exocyclic_oxygens
                    })

    num_monosaccharides = len(monosaccharide_rings)

    if num_monosaccharides != 4:
        return False, f"Found {num_monosaccharides} monosaccharide units, need exactly 4"

    # Build connectivity between monosaccharide rings via glycosidic bonds
    ring_connections = {}
    for idx_a, ring_a in enumerate(monosaccharide_rings):
        ring_a_atoms = ring_a['ring_atoms']
        exocyclic_oxygens = ring_a['exocyclic_oxygens']
        for (ring_carbon_idx, oxygen_idx) in exocyclic_oxygens:
            oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
            # Neighbors of oxygen atom
            neighbors = oxygen_atom.GetNeighbors()
            for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == ring_carbon_idx:
                    continue  # Skip the ring carbon it's attached to
                if neighbor.GetAtomicNum() == 6:
                    # Check if neighbor atom is in another monosaccharide ring
                    for idx_b, ring_b in enumerate(monosaccharide_rings):
                        if idx_b == idx_a:
                            continue
                        if neighbor_idx in ring_b['ring_atoms']:
                            # Found a glycosidic bond between ring_a and ring_b
                            # Record the connection
                            ring_connections.setdefault(idx_a, set()).add(idx_b)
                            ring_connections.setdefault(idx_b, set()).add(idx_a)

    # Check connectivity of monosaccharide rings
    def find_connected_components(connections):
        visited = set()
        components = []

        for node in connections:
            if node not in visited:
                stack = [node]
                component = set()
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        component.add(current)
                        neighbors = connections.get(current, [])
                        for neighbor in neighbors:
                            if neighbor not in visited:
                                stack.append(neighbor)
                components.append(component)
        return components

    components = find_connected_components(ring_connections)
    # Look for a connected component with exactly 4 monosaccharide rings
    for component in components:
        if len(component) == 4:
            return True, "Contains four monosaccharide units linked via glycosidic bonds"

    return False, "Monosaccharide units are not properly connected via glycosidic bonds"
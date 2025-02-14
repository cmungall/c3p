"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3β-hydroxy steroid whose skeleton is closely related to cholestan-3-ol,
    allowing for additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()
    number_of_rings = len(atom_rings)

    if number_of_rings < 4:
        return False, "Molecule has fewer than 4 rings"

    # Build list of rings with their sizes and atom indices
    rings = []
    for idx, atom_ring in enumerate(atom_rings):
        rings.append({'index': idx, 'size': len(atom_ring), 'atoms': set(atom_ring)})

    # Build ring adjacency (fusion) information
    ring_adjacency = {i: [] for i in range(number_of_rings)}
    for i in range(number_of_rings):
        for j in range(i+1, number_of_rings):
            # Rings are fused if they share two or more atoms (a bond)
            shared_atoms = rings[i]['atoms'] & rings[j]['atoms']
            if len(shared_atoms) >= 2:
                ring_adjacency[i].append(j)
                ring_adjacency[j].append(i)

    # Search for sterol skeleton: fused rings of sizes 6-6-6-5
    steroid_found = False
    for i in range(number_of_rings):
        if rings[i]['size'] != 6:
            continue
        for j in ring_adjacency[i]:
            if rings[j]['size'] != 6:
                continue
            for k in ring_adjacency[j]:
                if rings[k]['size'] != 6:
                    continue
                if k == i or k not in ring_adjacency[j]:
                    continue
                for l in ring_adjacency[k]:
                    if rings[l]['size'] != 5:
                        continue
                    if l == i or l == j or l == k or l not in ring_adjacency[k]:
                        continue
                    # Confirm that rings are fused in order i-j, j-k, k-l
                    if j in ring_adjacency[i] and k in ring_adjacency[j] and l in ring_adjacency[k]:
                        steroid_found = True
                        ringA_atoms = rings[i]['atoms']
                        ringB_atoms = rings[j]['atoms']
                        ringC_atoms = rings[k]['atoms']
                        ringD_atoms = rings[l]['atoms']
                        break
                if steroid_found:
                    break
            if steroid_found:
                break
        if steroid_found:
            break

    if not steroid_found:
        return False, "No steroid backbone found (cyclopentanoperhydrophenanthrene)"

    # Check for 3β-hydroxyl group on ring A
    # Assuming ring A is the first six-membered ring (rings[i])
    has_3beta_hydroxyl = False
    for atom_idx in ringA_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    # Found an oxygen atom (hydroxyl group) attached to carbon
                    has_3beta_hydroxyl = True
                    break
            if has_3beta_hydroxyl:
                break

    if not has_3beta_hydroxyl:
        return False, "No 3β-hydroxyl group found on the steroid backbone"

    return True, "Contains steroid backbone with 3β-hydroxyl group"
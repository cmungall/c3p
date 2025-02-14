"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:17754 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo().AtomRings()

    # Collect candidate sugar rings (5 or 6 membered rings with one oxygen, not connected to aromatic systems)
    candidate_rings = []
    for ring in rings:
        # Check ring size (5 or 6 members)
        if len(ring) == 5 or len(ring) == 6:
            # Count oxygen atoms in the ring
            num_oxygen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if num_oxygen == 1:
                # Ensure other atoms are carbons
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() in [6, 8] for idx in ring):
                    # Exclude rings connected to aromatic systems
                    ring_has_aromatic_connection = False
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIdx() not in ring and neighbor.GetIsAromatic():
                                ring_has_aromatic_connection = True
                                break
                        if ring_has_aromatic_connection:
                            break
                    if ring_has_aromatic_connection:
                        continue  # Skip this ring
                    candidate_rings.append(ring)

    amino_group_found = False  # Initialize here

    if candidate_rings:
        # For each candidate ring, check for amino groups replacing hydroxy groups
        for ring in candidate_rings:
            amino_group_found = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    # Check if carbon is sp3 hybridized
                    if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                        continue
                    # Check neighbors for nitrogen atoms
                    neighbors = atom.GetNeighbors()
                    for neighbor in neighbors:
                        if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                            # Check if bond is single (to exclude amides, etc.)
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                amino_group_found = True
                                break
                    if amino_group_found:
                        break
            if amino_group_found:
                return True, "Contains sugar ring with amino group(s) replacing hydroxy group(s)"

    # If no candidate rings with amino groups, proceed to check for acyclic amino sugars

    # Now, check for acyclic amino sugar patterns
    # Find all carbons with OH or NH2 groups
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    carbons_with_oh_or_nh2 = []
    for carbon in carbons:
        has_oh_or_nh2 = False
        has_amino_group = False
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if O is hydroxyl group
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetTotalDegree() == 1:
                    has_oh_or_nh2 = True
            elif neighbor.GetAtomicNum() == 7:
                # Check if N is amino group
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetTotalDegree() <= 3:
                    has_oh_or_nh2 = True
                    has_amino_group = True
        if has_oh_or_nh2:
            carbons_with_oh_or_nh2.append((carbon, has_amino_group))

    # Look for chains of such carbons
    from collections import deque
    visited = set()
    for carbon, has_amino in carbons_with_oh_or_nh2:
        if carbon.GetIdx() in visited:
            continue
        queue = deque()
        queue.append((carbon, has_amino, 1))
        visited.add(carbon.GetIdx())
        amino_in_chain = has_amino
        while queue:
            current_atom, current_has_amino, chain_length = queue.popleft()
            if chain_length >= 4 and amino_in_chain:
                return True, "Contains acyclic amino sugar"
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    if neighbor.GetIdx() in visited:
                        continue
                    # Check if neighbor has OH or NH2 group
                    neighbor_has_oh_or_nh2 = False
                    neighbor_has_amino = False
                    for nb in neighbor.GetNeighbors():
                        if nb.GetIdx() == current_atom.GetIdx():
                            continue
                        if nb.GetAtomicNum() == 8:
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nb.GetIdx())
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nb.GetTotalDegree() == 1:
                                neighbor_has_oh_or_nh2 = True
                        elif nb.GetAtomicNum() == 7:
                            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nb.GetIdx())
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nb.GetTotalDegree() <=3:
                                neighbor_has_oh_or_nh2 = True
                                neighbor_has_amino = True
                    if neighbor_has_oh_or_nh2:
                        visited.add(neighbor.GetIdx())
                        queue.append((neighbor, neighbor_has_amino, chain_length+1))
                        amino_in_chain = amino_in_chain or neighbor_has_amino

    return False, "Does not contain amino sugar ring or acyclic amino sugar"
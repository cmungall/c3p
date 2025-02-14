"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is an amino acid whose side chain is capable of forming one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify alpha carbon(s)
    alpha_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            n_neighbor = False
            carboxyl_neighbor = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 7:
                    n_neighbor = True
                elif nbr.GetAtomicNum() == 6:
                    # Potential carboxyl carbon
                    nbr_c = nbr
                    o_count = 0
                    for nbr_c_nbr in nbr_c.GetNeighbors():
                        if nbr_c_nbr.GetAtomicNum() == 8:
                            o_count += 1
                    if o_count >= 2:
                        carboxyl_neighbor = True
            if n_neighbor and carboxyl_neighbor:
                alpha_carbons.append(atom)

    if not alpha_carbons:
        return False, "No alpha carbon connected to both amino and carboxyl groups found"

    # Assume first alpha carbon
    alpha_carbon = alpha_carbons[0]
    alpha_carbon_idx = alpha_carbon.GetIdx()

    # Identify backbone atoms
    backbone_indices = set()
    backbone_indices.add(alpha_carbon_idx)

    # Identify nitrogen neighbor
    nitrogen_idx = None
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 7:
            nitrogen_idx = nbr.GetIdx()
            backbone_indices.add(nitrogen_idx)
            break
    if nitrogen_idx is None:
        return False, "Alpha carbon does not have a nitrogen neighbor"

    # Identify carboxyl carbon neighbor
    carboxyl_c_idx = None
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != nitrogen_idx:
            # Check if this carbon is connected to two oxygens
            nbr_c = nbr
            oxygen_indices = []
            for nbr_c_nbr in nbr_c.GetNeighbors():
                if nbr_c_nbr.GetAtomicNum() == 8:
                    oxygen_indices.append(nbr_c_nbr.GetIdx())
            if len(oxygen_indices) >=2:
                carboxyl_c_idx = nbr_c.GetIdx()
                backbone_indices.add(carboxyl_c_idx)
                backbone_indices.update(oxygen_indices)
                break

    if carboxyl_c_idx is None:
        return False, "Alpha carbon does not have a carboxyl group neighbor"

    # Identify side chain start atoms
    side_chain_start_atoms = []
    for nbr in alpha_carbon.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx not in backbone_indices:
            side_chain_start_atoms.append(nbr_idx)
            
    if not side_chain_start_atoms:
        return False, "No side chain found"

    # Get side chain atoms
    from collections import deque

    def get_side_chain_atoms(mol, side_chain_start_atoms, backbone_indices):
        side_chain_atoms = set()
        visited = set()
        queue = deque(side_chain_start_atoms)
        while queue:
            atom_idx = queue.popleft()
            if atom_idx in backbone_indices or atom_idx in visited:
                continue
            visited.add(atom_idx)
            side_chain_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in backbone_indices and nbr_idx not in visited:
                    queue.append(nbr_idx)
        return side_chain_atoms

    side_chain_atoms = get_side_chain_atoms(mol, side_chain_start_atoms, backbone_indices)

    # Check for hydrogen bond donors and acceptors in side chain
    hbond_donor_atoms = []
    hbond_acceptor_atoms = []

    for idx in side_chain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        atomic_num = atom.GetAtomicNum()
        # Check for hydrogen bond donors
        if atomic_num in [7, 8, 16]:  # N, O, S
            num_H = 0
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 1:  # Hydrogen
                    num_H += 1
            if num_H > 0:
                hbond_donor_atoms.append(idx)
        # Check for hydrogen bond acceptors
        if atomic_num in [7, 8, 16]:  # N, O, S
            # Exclude positively charged atoms (e.g., protonated amines)
            if atom.GetFormalCharge() >= 0:
                hbond_acceptor_atoms.append(idx)

    if hbond_donor_atoms or hbond_acceptor_atoms:
        return True, "Side chain can form hydrogen bonds"
    else:
        return False, "Side chain cannot form hydrogen bonds"
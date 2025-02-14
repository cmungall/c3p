"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: hydroxy fatty acid
"""

from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is any fatty acid carrying one or more hydroxy substituents.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule contains only C, H, and O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, f"Contains heteroatoms other than C, H, and O: {atom.GetSymbol()}"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, which is not typical for fatty acids"
    
    # Check for carboxylic acid group (both protonated and deprotonated forms)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O-,OH]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"
    
    # Get carboxylic acid carbon indices
    carboxy_carbons = [match[0] for match in carboxy_matches]
    
    # Check for aliphatic chain connected to carboxylic acid carbon
    max_chain_length = 0
    longest_chain = []
    for c_idx in carboxy_carbons:
        visited = set()
        chain = []
        def dfs(atom_idx, length):
            nonlocal max_chain_length, longest_chain
            if atom_idx in visited:
                return
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            # Skip carboxylic acid oxygens
            if atom.GetAtomicNum() != 6:
                return
            chain.append(atom_idx)
            if length > max_chain_length:
                max_chain_length = length
                longest_chain = chain.copy()
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
                if not bond.IsInRing():
                    dfs(n_idx, length + 1)
            chain.pop()
            visited.remove(atom_idx)
        for neighbor in mol.GetAtomWithIdx(c_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Only proceed with carbon atoms
                dfs(neighbor.GetIdx(), 1)
    if max_chain_length < 4:
        return False, f"Aliphatic chain length is {max_chain_length}, which is too short for a fatty acid"
    
    # Check for hydroxy groups (excluding carboxylic acid hydroxyl)
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    # Get indices of hydroxy oxygens
    hydroxy_oxygens = [match[0] for match in hydroxy_matches]
    # Remove carboxylic acid oxygens from hydroxy list
    carboxy_oxygens = [match[2] for match in carboxy_matches]
    hydroxy_oxygens = [idx for idx in hydroxy_oxygens if idx not in carboxy_oxygens]
    if not hydroxy_oxygens:
        return False, "No hydroxy substituents found excluding carboxylic acid group"
    
    # Check if hydroxy groups are attached to the aliphatic chain
    hydroxy_on_chain = False
    chain_atom_set = set(longest_chain)
    for o_idx in hydroxy_oxygens:
        oxygen_atom = mol.GetAtomWithIdx(o_idx)
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in chain_atom_set:
                hydroxy_on_chain = True
                break
        if hydroxy_on_chain:
            break
    if not hydroxy_on_chain:
        return False, "Hydroxy groups are not on the aliphatic chain"
    
    return True, "Molecule is a hydroxy fatty acid"
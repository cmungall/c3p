"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def get_carbon_chain(mol, start_idx, end_idx=None):
    """Find a continuous carbon chain between start and end atoms"""
    if end_idx is None:
        # Find the furthest carbon atom
        visited = set()
        queue = deque([(start_idx, [start_idx])])
        longest_path = []
        
        while queue:
            current, path = queue.popleft()
            if len(path) > len(longest_path):
                longest_path = path
                
            atom = mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and 
                    neighbor.GetIdx() not in visited and 
                    neighbor.GetIdx() not in path):
                    new_path = path + [neighbor.GetIdx()]
                    queue.append((neighbor.GetIdx(), new_path))
                    visited.add(neighbor.GetIdx())
        
        return longest_path
    else:
        # Find shortest path between start and end
        return Chem.GetShortestPath(mol, start_idx, end_idx)

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[OX2H][CX3](=[OX1])")
    cooh_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if not cooh_matches:
        return False, "No carboxylic acid group found"
    if len(cooh_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Get the carboxylic acid carbon
    cooh_carbon_idx = cooh_matches[0][1]
    
    # Find terminal hydroxyl groups
    terminal_oh_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    
    if not oh_matches:
        return False, "No hydroxyl groups found"
    
    # Find valid terminal hydroxyl groups (with only one carbon neighbor)
    terminal_carbons = []
    for match in oh_matches:
        carbon_idx = match[0]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_carbons.append(carbon_idx)
    
    if not terminal_carbons:
        return False, "No terminal hydroxyl groups found"
    
    # Check ring characteristics
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
        if len(ring) > 3:  # Allow only epoxy rings
            return False, "Contains large ring structure"
    
    # For each terminal OH, check if it forms a valid fatty acid chain
    valid_chain_found = False
    for terminal_carbon in terminal_carbons:
        chain = get_carbon_chain(mol, cooh_carbon_idx, terminal_carbon)
        
        if chain and len(chain) >= 6:  # Minimum chain length for fatty acids
            # Count hydroxyl groups along the chain
            oh_count = 0
            for idx in chain:
                atom = mol.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 2:
                        oh_count += 1
            
            # Check OH density (exclude sugar-like structures)
            if oh_count / len(chain) > 0.5:
                continue
            
            # Check for excessive branching
            has_excessive_branching = False
            for idx in chain:
                if idx in ring_atoms:
                    continue
                atom = mol.GetAtomWithIdx(idx)
                carbon_neighbors = len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])
                if carbon_neighbors > 2:
                    has_excessive_branching = True
                    break
            
            if not has_excessive_branching:
                valid_chain_found = True
                break
    
    if not valid_chain_found:
        return False, "No valid fatty acid chain found with terminal hydroxyl"
    
    # Additional checks
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains elements other than C, H, O"
    
    return True, "Contains omega-hydroxyl fatty acid structure with appropriate chain length"
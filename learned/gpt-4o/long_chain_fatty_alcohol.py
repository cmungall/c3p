"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a primary alcohol group and a carbon chain length of 13 to 22 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary alcohol group (OH); must be bound to a carbon atom
    primary_alcohol_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol group (C-OH) found"
    
    # Traverse the molecule to find the longest carbon chain
    atom_types = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if not any(atom == 6 for atom in atom_types):  # ensure there are carbon atoms
        return False, "No carbon atoms present"

    # Initialize the longest carbon chain length
    longest_chain_length = 0
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Only start from carbon atoms
            # Explore the depth of carbon chains using BFS
            visited = {atom.GetIdx()}
            queue = [(atom, 0)]
            max_depth = 0
            
            while queue:
                current_atom, depth = queue.pop(0)
                max_depth = max(max_depth, depth)
                
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                        visited.add(neighbor.GetIdx())
                        queue.append((neighbor, depth + 1))

            # If the longest chain in this path is greater than previously found, update it
            if max_depth + 1 > longest_chain_length:
                longest_chain_length = max_depth + 1

    # Validate chain length
    if 13 <= longest_chain_length <= 22:
        return True, f"Contains a primary alcohol group and a suitable carbon chain length of {longest_chain_length}."

    return False, f"Longest carbon chain length is {longest_chain_length}, outside the 13-22 range"
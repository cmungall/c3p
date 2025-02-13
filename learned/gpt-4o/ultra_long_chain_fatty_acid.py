"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is characterized by a carbon chain length greater than C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (COOH)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    def get_longest_linear_carbon_chain(mol):
        """ Finds the longest path of linear, acyclic carbon atoms in the molecule """
        longest_chain = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Only consider carbon atoms
                visited = set()
                current_chain = explore_linear_chain(atom, visited)
                longest_chain = max(longest_chain, current_chain)
        return longest_chain

    def explore_linear_chain(atom, visited):
        """ Explore acyclic carbon chain starting at the given atom """
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():  # Avoid aromatic carbons
            return 0
        visited.add(atom.GetIdx())
        max_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6 and len(neighbor.GetNeighbors()) <= 4:
                # Ensure acyclic by limiting bonds
                chain_length = explore_linear_chain(neighbor, visited.copy())
                max_length = max(max_length, chain_length)
        return max_length + 1

    # Calculate the longest linear carbon chain in the molecule
    longest_chain_length = get_longest_linear_carbon_chain(mol)
    
    if longest_chain_length > 27:
        return True, f"Contains a carbon chain of length {longest_chain_length}, qualifying as ultra-long-chain"

    return False, f"Longest carbon chain is {longest_chain_length}; not greater than 27"
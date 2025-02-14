"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid which has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group (O=C-O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not matches:
        return False, "No carboxylic acid group found"

    # Use the first carboxylic acid group found
    carboxylic_c_idx = matches[0][0]  # Carbonyl carbon index

    # Identify terminal carbons (degree 1 carbons excluding the carboxylic carbon)
    terminal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            idx = atom.GetIdx()
            if idx == carboxylic_c_idx:
                continue  # Skip the carboxylic carbon
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            if len(neighbors) == 1:
                terminal_carbons.append(idx)

    # Initialize maximum chain length
    max_chain_length = 0

    # Compute the longest linear carbon chain starting from the carboxylic carbon
    for t_c_idx in terminal_carbons:
        # Get shortest path between carboxylic carbon and terminal carbon
        path = Chem.rdmolops.GetShortestPath(mol, carboxylic_c_idx, t_c_idx)
        carbons_in_path = [idx for idx in path if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

        # Check for branching in the chain
        branching = False
        for idx in carbons_in_path[1:-1]:  # Exclude end carbons
            atom = mol.GetAtomWithIdx(idx)
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            if len(neighbors) != 2:
                branching = True
                break
        if branching:
            continue  # Skip this path due to branching

        # Update maximum chain length
        chain_length = len(carbons_in_path)
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    if max_chain_length == 0:
        return False, "No linear carbon chain found"

    if max_chain_length > 22:
        return True, f"Chain length is {max_chain_length}, which is greater than 22"
    else:
        return False, f"Chain length is {max_chain_length}, which is not greater than 22"
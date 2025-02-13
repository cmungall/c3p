"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a carbon chain length greater than 22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    max_chain_length = 0
    
    def count_longest_chain_from_atom(atom_index):
        visited = set()
        stack = [(atom_index, 0)]  # (current_atom_index, current_chain_length)
        max_length = 0

        while stack:
            current_atom, current_length = stack.pop()
            if current_atom not in visited:
                visited.add(current_atom)
                found_extension = False
                
                for neighbor in mol.GetAtomWithIdx(current_atom).GetNeighbors():
                    neighbor_index = neighbor.GetIdx()

                    # Continue if we find a carbon that is not the one in the carboxyl group and is not visited
                    if neighbor_index not in visited and neighbor.GetAtomicNum() == 6:  # is a carbon
                        stack.append((neighbor_index, current_length + 1))
                        found_extension = True

                if not found_extension:
                    max_length = max(max_length, current_length)
        
        return max_length

    # Measure max chain length from each carboxylic acid group's primary carbon atom
    for match in carboxyl_matches:
        carbon_in_carboxyl = match[0]  # The first atom in the match is the carbon atom
        chain_length = count_longest_chain_from_atom(carbon_in_carboxyl)

        if chain_length > max_chain_length:
            max_chain_length = chain_length

    # Check chain length criteria for very long-chain fatty acid (greater than 22 carbons)
    if max_chain_length > 22:
        return True, f"Longest carbon chain length is {max_chain_length}, which is greater than 22"

    return False, f"Longest carbon chain length is {max_chain_length}, not greater than 22"
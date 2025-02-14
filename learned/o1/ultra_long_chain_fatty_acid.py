"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: ultra-long-chain fatty acid
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is defined as any very long-chain fatty acid which has a chain length greater than C27.

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

    # Look for exactly one carboxylic acid functional group [CX3](=O)[O;H1,-1]
    # This pattern matches both protonated and deprotonated carboxylic acids
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) != 1:
        return False, f"Expected one carboxylic acid group, found {len(matches)}"

    # Check for rings larger than cyclopropane (size > 3)
    ring_info = mol.GetRingInfo()
    has_large_ring = False
    for ring in ring_info.AtomRings():
        if len(ring) > 3:
            has_large_ring = True
            break
    if has_large_ring:
        return False, "Molecule contains rings larger than cyclopropane"

    # Find the carboxyl carbon index
    carboxyl_carbon_idx = matches[0][0]

    # Function to find the longest carbon chain connected to a given atom
    def find_longest_chain(atom_idx, visited):
        visited.add(atom_idx)
        max_length = 1  # Include current atom
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            # Only consider carbon atoms not visited yet
            if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                length = 1 + find_longest_chain(neighbor_idx, visited)
                if length > max_length:
                    max_length = length
        return max_length

    visited = set()
    chain_length = find_longest_chain(carboxyl_carbon_idx, visited)

    if chain_length > 27:
        return True, f"Longest carbon chain is {chain_length} carbons"
    else:
        return False, f"Longest carbon chain is {chain_length} carbons, which is not greater than 27"
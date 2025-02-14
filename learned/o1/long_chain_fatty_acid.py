"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:15904 long-chain fatty acid
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as a fatty acid with a chain length ranging from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a fatty acid"

    # Define SMARTS patterns for carboxylic acid and carboxylate groups
    carboxylic_acid_smarts = '[CX3](=O)[OX1H0-,OX2H1]'
    carboxylic_acid = Chem.MolFromSmarts(carboxylic_acid_smarts)
    matches = mol.GetSubstructMatches(carboxylic_acid)

    # Check for exactly one carboxylic acid group
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, expected exactly 1"

    # Get the carboxyl carbon index
    carboxyl_carbon_idx = matches[0][0]

    # Get the chain length from the carboxyl carbon
    chain_length = get_linear_chain_length(mol, carboxyl_carbon_idx)

    # Check if the chain length is within the specified range
    if 13 <= chain_length <= 22:
        return True, f"Chain length is {chain_length}, within 13-22"
    else:
        return False, f"Chain length is {chain_length}, not within 13-22"

def get_linear_chain_length(mol, carboxyl_carbon_idx):
    """
    Finds the length of the unbranched carbon chain starting from the carboxyl carbon.

    Args:
        mol (Mol): RDKit molecule object
        carboxyl_carbon_idx (int): Atom index of the carboxyl carbon

    Returns:
        int: Length of the carbon chain including the carboxyl carbon
    """
    chain_length = 1  # Include carboxyl carbon
    visited = set()
    current_atom_idx = carboxyl_carbon_idx
    previous_atom_idx = None

    while True:
        visited.add(current_atom_idx)
        current_atom = mol.GetAtomWithIdx(current_atom_idx)

        # Get neighboring carbon atoms excluding the previous atom
        neighbors = [
            neighbor for neighbor in current_atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != previous_atom_idx
        ]

        # Exclude carbons already visited
        neighbors = [neighbor for neighbor in neighbors if neighbor.GetIdx() not in visited]

        # If more than one neighbor (excluding previous), branching occurs
        if len(neighbors) > 1:
            break  # Branching occurs, so chain is no longer linear

        if len(neighbors) == 0:
            break  # End of chain

        # Move to the next carbon
        next_atom = neighbors[0]
        chain_length += 1
        previous_atom_idx = current_atom_idx
        current_atom_idx = next_atom.GetIdx()

    return chain_length
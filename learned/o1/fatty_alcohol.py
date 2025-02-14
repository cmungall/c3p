"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify hydroxyl groups attached to sp3 carbons (exclude phenols)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;!$(C[O,N,S,P]=O)]-[OX2H]")
    matches = mol.GetSubstructMatches(hydroxyl_pattern)

    if not matches:
        return False, "No suitable hydroxyl groups attached to sp3 carbons"

    # Iterate over each hydroxyl group
    for match in matches:
        carbon_idx = match[0]
        oxygen_idx = match[1]

        # Check if the carbon is in a ring
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if carbon_atom.IsInRing():
            continue  # Skip if carbon is part of a ring

        # Perform BFS to find all connected carbons (excluding rings)
        visited = set()
        queue = [carbon_idx]
        carbon_count = 0

        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)

            # Count only carbon atoms
            if atom.GetAtomicNum() == 6:
                carbon_count += 1
            else:
                continue  # Skip non-carbon atoms

            # Add neighboring carbons that are not in rings
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and not neighbor.IsInRing():
                    if neighbor.GetAtomicNum() == 6:
                        queue.append(neighbor_idx)

        if carbon_count >= 3:
            return True, f"Contains aliphatic alcohol with a chain of {carbon_count} carbons"

    return False, "No aliphatic chain with sufficient length connected to hydroxyl group"

# Example usage:
# smiles_list = [
#     "CCCCCCCCCCCCCC(O)CCC",  # tetradecan-5-ol
#     "CCCCCCCCCCCCCO",        # dodecan-1-ol
#     "CC(C)(O)C(=O)Cc1cnc2ccccc12",  # Should be False
# ]
# for smiles in smiles_list:
#     result, reason = is_fatty_alcohol(smiles)
#     print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")
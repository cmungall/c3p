"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has a terminal aldehyde group connected to a longer carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldehyde group pattern: R-CH=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group found"

    # Check terminal positioning and chain length for each aldehyde group match
    for match in matches:
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])

        # Count carbons in the connected substructure starting from the aldehyde carbon
        neighbors = [nbr for nbr in aldehyde_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(neighbors) != 1:
            continue
        
        chain_carbon = neighbors[0]
        visited = {aldehyde_carbon.GetIdx()}
        carbon_count = 0
        stack = [chain_carbon]

        while stack:
            current = stack.pop()
            if current.GetIdx() in visited or current.IsInRing():
                continue
            visited.add(current.GetIdx())
            if current.GetAtomicNum() == 6:
                carbon_count += 1
                for nbr in current.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                        stack.append(nbr)

        # A fatty aldehyde typically has a terminal chain longer than 5 carbons
        if carbon_count >= 5:
            return True, f"Has terminal aldehyde group in acyclic carbon chain of {carbon_count} carbons"

    return False, "Aldehyde group not terminal in proper chain structure"
"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is a fatty acyl-CoA where the fatty acyl group has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acyl thioester pattern (C(=O)S)
    acyl_thioester_pattern = Chem.MolFromSmarts('[#6](=O)S')

    # Find the acyl thioester group
    matches = mol.GetSubstructMatches(acyl_thioester_pattern)
    if not matches:
        return False, "No acyl-CoA thioester group found"

    # Assume the first match is the acyl thioester bond
    match = matches[0]
    carbonyl_carbon_idx = match[0]
    sulfur_idx = match[1]

    # Get the carbon attached to the carbonyl carbon (start of the fatty acyl chain)
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
    neighbors = [atom for atom in carbonyl_carbon.GetNeighbors() if atom.GetAtomicNum() not in (8, 16)]  # Exclude oxygen and sulfur
    if not neighbors:
        return False, "No fatty acyl chain attached to the carbonyl carbon"
    acyl_chain_start_atom = neighbors[0]
    acyl_chain_start_idx = acyl_chain_start_atom.GetIdx()

    # Traverse the fatty acyl chain and count carbons
    visited = set([carbonyl_carbon_idx, sulfur_idx])
    to_visit = [acyl_chain_start_idx]
    carbon_count = 1  # Start with 1 for the acyl_chain_start_atom

    while to_visit:
        current_idx = to_visit.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        if current_atom.GetAtomicNum() == 6:
            # Count the carbon
            carbon_count += 1
            # Add neighbors that are carbons or oxygens (allowing for double bonds and hydroxy groups)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() in (6, 8):
                    # Exclude atoms that would lead back to CoA (e.g., nitrogen, phosphorus)
                    to_visit.append(neighbor_idx)
        else:
            # Stop traversal if we hit heteroatoms other than oxygen (e.g., nitrogen in CoA)
            continue

    # Check if the chain length is greater than 22 carbons
    if carbon_count > 22:
        return True, f"Fatty acyl chain length is {carbon_count} carbons, which is greater than 22"
    else:
        return False, f"Fatty acyl chain length is {carbon_count} carbons, which is not greater than 22"
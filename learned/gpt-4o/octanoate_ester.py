"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is characterized by the ester group (-C(=O)O-) where the acid part 
    is specifically octanoic acid (8 carbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester group pattern (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found"
    
    # For each ester match, check for 8-carbon chain attached to the carbonyl
    for match in ester_matches:
        # Identify the carbonyl carbon
        carbonyl_c = match[0]

        # Check if there is a linear chain of exactly 8 carbons stemming from the carbonyl carbon
        chain_carbons = set([carbonyl_c])  # Initialize with carbonyl carbon
        queue = [carbonyl_c]

        while queue and len(chain_carbons) <= 8:
            current_carbon = queue.pop(0)

            for neighbor in mol.GetAtomWithIdx(current_carbon).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                
                if neighbor_idx not in chain_carbons and neighbor.GetAtomicNum() == 6:
                    chain_carbons.add(neighbor_idx)
                    queue.append(neighbor_idx)

        # Check if the chain has exactly 8 carbons (including carbonyl carbon)
        if len(chain_carbons) == 8:
            return True, "Contains octanoate ester group"
    
    return False, "Ester group found, but no octanoate (8-carbon chain)"
"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid with a carbon chain greater than C22.
    Ultra-long-chain fatty acids are a subset with chain length > C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group or its deprotonated form
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    

    max_chain_length = 0
    for acid_match in acid_matches:
        carbonyl_carbon_idx = acid_match[0]
        
        # Get the carbon adjacent to the carbonyl carbon
        neighbors = mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors()
        chain_start_carbons = [neighbor for neighbor in neighbors if neighbor.GetSymbol() == 'C']
        
        current_max_len = 0
        for start_carbon in chain_start_carbons:
            visited_atoms = {carbonyl_carbon_idx, start_carbon.GetIdx()} # keep track of visited atoms including carbonyl
            current_chain_length = 1
            
            
            queue = [(start_carbon.GetIdx(), 1, visited_atoms)] # start from the neighboring C, keeping track of visited atoms, and current length
           
            
            while queue:
              current_idx, current_len, current_visited = queue.pop(0)
              
              max_len_from_here = current_len
              
              
              neighbors = mol.GetAtomWithIdx(current_idx).GetNeighbors()
              
              next_carbons = []
              for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                if neighbor.GetSymbol() == "C" and neighbor_idx not in current_visited:
                      
                  ring_info = mol.GetRingInfo()
                  if not ring_info.IsAtomInRingOfSize(neighbor_idx, 3) and not ring_info.IsAtomInRingOfSize(neighbor_idx, 4) and not ring_info.IsAtomInRingOfSize(neighbor_idx, 5) and not ring_info.IsAtomInRingOfSize(neighbor_idx, 6):
                      next_carbons.append( (neighbor_idx, current_len + 1, set(current_visited) | {neighbor_idx}) )
              
              
              
              if not next_carbons: #we reached the end of the chain
                
                current_max_len = max(current_max_len, current_len)
              else:
                  queue.extend(next_carbons)
            
            current_max_len += 1 # account for carbonyl C
            max_chain_length = max(max_chain_length, current_max_len)



    if max_chain_length > 27:
        return True, f"Ultra-long-chain fatty acid (C={max_chain_length})"
    elif max_chain_length > 22:
        return True, f"Very long-chain fatty acid (C={max_chain_length})"
    else:
        return False, f"Not a very long-chain fatty acid (C={max_chain_length})"
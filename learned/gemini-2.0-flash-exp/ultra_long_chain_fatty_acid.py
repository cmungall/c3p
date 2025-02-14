"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid
    (chain length > C27) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a single carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(matches) != 1:
         return False, "Molecule does not contain exactly one carboxylic acid group."
    
    # Get the carbon atom index of the carboxyl carbon
    carboxyl_carbon_idx = matches[0][0]

    # Function to find the longest chain of carbons starting from carboxyl carbon
    def find_longest_chain(start_atom_idx, visited_atoms):
      
      max_len = 0
      for neighbor in mol.GetAtomWithIdx(start_atom_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_atoms:
          max_len = max(max_len, 1 + find_longest_chain(neighbor.GetIdx(), visited_atoms | {neighbor.GetIdx()}))

      return max_len

    longest_chain = find_longest_chain(carboxyl_carbon_idx, {carboxyl_carbon_idx}) +1 # +1 for the starting carbon

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)


    # Filter out non-fatty acid molecules using oxygen counts
    if o_count > 5:
        return False, f"Too many oxygen atoms for a fatty acid ({o_count})"

    if longest_chain > 27:
        return True, f"Longest carbon chain has {longest_chain} carbons, which is > 27."
    else:
        return False, f"Longest carbon chain has {longest_chain} carbons, which is <= 27."
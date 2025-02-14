"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA derived from a branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for CoA
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
          return False, "CoA moiety not found"

    # Define SMARTS pattern for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, "Incorrect number of thioester linkages"

    # Find carbonyl carbon of thioester
    carbonyl_carbon_indices = [match[0] for match in mol.GetSubstructMatches(thioester_pattern)]
    if not carbonyl_carbon_indices:
       return False, "Carbonyl group of thioester not found"
    
    carbonyl_carbon_index = carbonyl_carbon_indices[0]
    
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_index)
    
    # Get neighbors of carbonyl carbon
    neighbors = [atom.GetIdx() for atom in carbonyl_carbon.GetNeighbors()]
    
    # Find the carbon that is part of the fatty acyl chain
    fatty_chain_carbon_index = None
    for neighbor_index in neighbors:
      neighbor = mol.GetAtomWithIdx(neighbor_index)
      if neighbor.GetSymbol() == 'C':
        fatty_chain_carbon_index = neighbor_index
        break
    
    if fatty_chain_carbon_index is None:
      return False, "Fatty chain carbon not found"

    
    # Perform depth-first search, keep track of visited atoms and chain atoms
    visited_atoms = set()
    acyl_chain_atoms = set()

    def dfs(atom_index):
        if atom_index in visited_atoms:
            return False # avoid loops
        
        visited_atoms.add(atom_index)

        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() != 'C':
          return False
        
        acyl_chain_atoms.add(atom_index)
        
        
        
        carbon_neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'C']
        
        #Check for branching only on the chain. Only if > 2 carbon neighbors in acyl chain
        branching_neighbors = 0
        for neighbor_index in carbon_neighbors:
           if neighbor_index in acyl_chain_atoms or neighbor_index not in visited_atoms:
                branching_neighbors+=1
        if branching_neighbors > 2:
            return True  # Found branch on acyl chain
        

        for neighbor_index in carbon_neighbors:
          if dfs(neighbor_index):
            return True
        
        return False  # no branch found yet

    is_branched = dfs(fatty_chain_carbon_index)
    
    if is_branched:
        return True, "Molecule is a branched-chain fatty acyl-CoA"
    else:
      return False, "No branch found on the fatty acyl chain."
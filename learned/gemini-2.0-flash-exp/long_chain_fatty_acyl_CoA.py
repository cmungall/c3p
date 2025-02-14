"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a coenzyme A derivative with a fatty acid chain of 13-22 carbons,
    linked via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None # Invalid SMILES

    # Check for CoA substructure (adenosine diphosphate part + pantetheine fragment)
    # Simplified SMARTS pattern for CoA part, focusing on the phosphate and nucleotide region
    coa_pattern = Chem.MolFromSmarts("COP(=O)([OX1])OP(=O)([OX1])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Molecule does not contain the required coenzyme A substructure"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=[OX1])")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
         return False, "Molecule does not contain a thioester group"
    
    # Find the carbon attached to the carbonyl of the thioester
    thioester_carbon = None
    for match in thioester_matches:
      for atom_index in match:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2: #Carbonyl carbon
           thioester_carbon = atom
           break
      if thioester_carbon:
        break
    
    if thioester_carbon is None:
      return False, "Could not identify the fatty acid chain"
    
    # Get the atoms connected to the carbonyl carbon, skipping the sulfur
    neighbor_atoms = [neighbor for neighbor in thioester_carbon.GetNeighbors() if neighbor.GetSymbol() != 'S']
    
    if len(neighbor_atoms) != 1:
        return False, "Incorrect number of connections to the carbonyl carbon"
    
    fatty_acid_start = neighbor_atoms[0]
    if fatty_acid_start.GetSymbol() != 'C':
        return False, "Fatty acid chain not directly linked to the thioester group"

    # Trace and count carbons in the fatty acid chain using DFS and exclude Coenzyme A atoms.
    visited_atoms = set()
    stack = [fatty_acid_start]
    carbon_count = 0
    
    coa_atoms = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        for idx in match:
            coa_atoms.add(idx)
            
    while stack:
      current_atom = stack.pop()
      if current_atom.GetIdx() in visited_atoms or current_atom.GetIdx() in coa_atoms:
        continue
      
      visited_atoms.add(current_atom.GetIdx())
      
      if current_atom.GetAtomicNum() == 6: # Carbon
        carbon_count+=1

      for neighbor in current_atom.GetNeighbors():
        stack.append(neighbor)
            
    if carbon_count < 13 or carbon_count > 22:
        return False, f"Fatty acid chain has {carbon_count} carbons, should be between 13 and 22."
    

    return True, "Molecule is a long-chain fatty acyl-CoA"
"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the core CoA structure (simplified to include key parts)
    # This is just a fragment and might not work well in all situations
    coa_pattern = Chem.MolFromSmarts("C[C@H](O)[C@@H](COP(=O)([O-])OP(=O)([O-])OCC1O[C@H]([C@H](O)[C@@H](O)1)n1cnc2c(N)ncnc12)C(C)(C)") 
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA substructure not found"
    
    # Check for thioester bond (-C(=O)-S-) connected to CoA substructure
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"
    
    # Find the acyl carbon chain attached to the thioester
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)[CX4]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "No acyl chain attached"

    # Count the carbons in acyl chain
    acyl_chain_start_atom_idx = acyl_chain_matches[0][1] # second atom in the match (the one directly attached to the carbonyl)
    
    # Get the heavy atoms attached to the acyl carbon, and traverse the chain
    chain_carbons = 0
    queue = [acyl_chain_start_atom_idx]
    visited = set()
    while queue:
       current = queue.pop(0)
       if current in visited:
           continue
       visited.add(current)
       atom = mol.GetAtomWithIdx(current)
       if atom.GetAtomicNum() == 6:
           chain_carbons+=1
           for neighbor in atom.GetNeighbors():
               if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                   queue.append(neighbor.GetIdx())
    
    # Remove the 1st acyl carbon that is part of the carbonyl
    chain_carbons -=1
    
    if  chain_carbons < 6 or chain_carbons > 12:
           return False, f"Acyl chain has {chain_carbons} carbons, must be between 6 and 12"

    # Verify the total charge of -4 using the formal charges of the atoms
    total_charge = 0
    for atom in mol.GetAtoms():
      total_charge += atom.GetFormalCharge()

    if total_charge != -4:
        return False, f"Molecule has {total_charge} charge instead of -4"

    return True, "Matches medium-chain fatty acyl-CoA(4-) criteria"
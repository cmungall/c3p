"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group 
    of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the Coenzyme A SMARTS pattern (simplified)
    coa_pattern = Chem.MolFromSmarts("NC(=O)CNC(=O)[C@H](O)C(C)(C)COP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A moiety found"
    
    # Define the thioester linkage SMARTS pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Find all carbon chains attached to the thiol group through thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    for match in thioester_matches:
        carbon_chain = []
        
        # Go through the bonded atoms to map the carbon chain
        atom_idx = match[0]  # The carbon of the thioester
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # Traverse the chain starting from the thioester carbon
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                chain_atom = neighbor
                break
        else:
            continue  # No carbon found, move to next thioester
        
        # Simple traversal to count carbons
        carbon_count = 1  # Including the initial carbon
        visited = set([atom_idx, match[1]])  # Start with thioester atoms
        to_visit = [chain_atom.GetIdx()]
        visited.add(chain_atom.GetIdx())
        
        while to_visit:
            current_atom_idx = to_visit.pop()
            carbon_count += 1
            current_atom = mol.GetAtomWithIdx(current_atom_idx)
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    to_visit.append(neighbor.GetIdx())
                    visited.add(neighbor.GetIdx())
        
        # Check if carbon chain length is within the desired range
        if 13 <= carbon_count <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {carbon_count} carbon atoms in fatty acid chain"
    
    return False, "No suitable long-chain fatty acid chain found"
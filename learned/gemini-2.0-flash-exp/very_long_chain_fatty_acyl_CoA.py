"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA (greater than C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove hydrogen atoms before sanitizing
    mol = Chem.RemoveHs(mol)
    
    #Sanitize to remove explicit charges
    Chem.SanitizeMol(mol)

    # Define the CoA substructure
    coa_substructure = Chem.MolFromSmiles('SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_substructure):
        return False, "CoA substructure not found"

    # Find the ester/thioester substructure (C(=O)S) to identify the carbonyl carbon connected to the fatty acid
    carbonyl_substructure = Chem.MolFromSmarts('[CX3](=O)[S]')
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_substructure)

    if len(carbonyl_matches) == 0:
       return False, "No carbonyl linked to sulfur found"
    
    if len(carbonyl_matches) > 1:
        return False, "Multiple carbonyls linked to sulfur found, cannot identify fatty acid chain"


    # Get the carbon atom index in the carbonyl group
    carbonyl_carbon_index = carbonyl_matches[0][0]

    # Get the indices of all atoms in the CoA substructure
    coa_atom_indices = set()
    for match in mol.GetSubstructMatches(coa_substructure):
        coa_atom_indices.update(match)
    
    # Find the directly attached carbon to the carbonyl
    neighboring_carbons = []
    for neighbor in mol.GetAtomWithIdx(carbonyl_carbon_index).GetNeighbors():
         if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in coa_atom_indices:
            neighboring_carbons.append(neighbor.GetIdx())

    if len(neighboring_carbons) == 0:
         return False, "No carbon found attached to the carbonyl group."
    
    if len(neighboring_carbons) > 1:
        return False, "Multiple carbons found attached to the carbonyl group, cannot identify fatty acid chain"
    
    first_chain_carbon_index = neighboring_carbons[0]
    
    # Now perform BFS from the first_chain_carbon_index to find the rest of the carbon chain.
    fatty_acid_carbon_count = 1  # Include the carbon attached to the carbonyl
    visited_atoms = {carbonyl_carbon_index, first_chain_carbon_index}
    queue = [first_chain_carbon_index]


    while queue:
        current_atom_index = queue.pop(0)
        
        # Count the number of directly attached carbons to the current atom, that are NOT part of the coA.
        # This is a proxy to make sure we are only counting the *chain* of carbons, and not branching.
        connected_carbons = 0
        next_carbon_index = None
        
        for neighbor in mol.GetAtomWithIdx(current_atom_index).GetNeighbors():
           neighbor_idx = neighbor.GetIdx()
           if neighbor.GetAtomicNum() == 6 and neighbor_idx not in coa_atom_indices and neighbor_idx not in visited_atoms:
                connected_carbons += 1
                next_carbon_index = neighbor_idx
        
        if connected_carbons == 1:
            fatty_acid_carbon_count += 1
            visited_atoms.add(next_carbon_index)
            queue.append(next_carbon_index)
        elif connected_carbons > 1:
            # Branched chain, stop counting
            break
    
    if fatty_acid_carbon_count > 22:
       return True, f"Fatty acid chain contains {fatty_acid_carbon_count} carbons, which is greater than 22"
    else:
       return False, f"Fatty acid chain contains {fatty_acid_carbon_count} carbons, which is not greater than 22"
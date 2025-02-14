"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is a fatty acid that contains one or more rings as part of the fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid or ester group
    acid_or_ester_pattern = Chem.MolFromSmarts('[CX3,CX4](=[OX1])[OX2]')
    acid_or_ester_matches = mol.GetSubstructMatches(acid_or_ester_pattern)

    if not acid_or_ester_matches:
        return False, "No carboxylic acid or ester group found"
    
    found_cyclic_chain = False
    for match in acid_or_ester_matches:
        for atom_index in match:
            if mol.GetAtomWithIdx(atom_index).GetAtomicNum() == 6: #carbon atom of acid/ester group
                
                chain_carbons = []
                
                #Function to trace carbon chain
                def find_chain(current_atom, visited_atoms):
                    
                    if current_atom.GetIdx() in visited_atoms:
                       return
                    
                    visited_atoms.add(current_atom.GetIdx())
                    
                    
                    
                    if current_atom.GetAtomicNum() == 6: #add this carbon to the chain
                       chain_carbons.append(current_atom.GetIdx())
                    
                    
                    for neighbor in current_atom.GetNeighbors():
                         if neighbor.GetAtomicNum() == 6: #Check if neighbour is also a carbon atom
                             find_chain(neighbor, visited_atoms)

                
                
                visited = set()
                
                
                for neighbor in mol.GetAtomWithIdx(atom_index).GetNeighbors(): #explore neighbours of the carbon atom of the acid group
                    if neighbor.GetAtomicNum() == 6:
                       find_chain(neighbor, visited) #find all carbons of the fatty acid chain.
                
                
                #remove carbons within rings to not count twice
                chain_carbons_no_ring = []
                for carbon_idx in chain_carbons:
                    carbon = mol.GetAtomWithIdx(carbon_idx)
                    in_ring = carbon.IsInRing()
                    if not in_ring:
                       chain_carbons_no_ring.append(carbon_idx)
                
                #Check length of fatty acid chain excluding rings.
                if len(chain_carbons_no_ring) < 4:
                   continue # continue to next acid group, this is not a fatty acid
            
                
                ring_info = mol.GetRingInfo()

                
                #check if a carbon of the chain is part of a ring
                for carbon_idx in chain_carbons:
                     carbon = mol.GetAtomWithIdx(carbon_idx)
                     if carbon.IsInRing():
                         found_cyclic_chain = True
                         break
                if found_cyclic_chain:
                     break #break the loop, we found a match
        if found_cyclic_chain:
              break #break the main loop


    if not found_cyclic_chain:
        return False, "No ring found as part of a fatty acid carbon chain"

    return True, "Contains ring and fatty acid substructures"
"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA consists of Coenzyme A linked to a fatty acid with 6-12 carbons via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for CoA core
    # The following SMARTS captures the key part of the coenzyme A moiety, including the pyrophosphate,
    # ribose, and adenine parts, while being flexible about the exact position of the thiol and other connecting parts.
    # Also, we include the amide that attaches to the pantetheine moiety.
    coa_pattern = Chem.MolFromSmarts('C[C@H](O)COP(=O)(OP(=O)(OCC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12)O)O.NC(=O)CCNC(=O)CCS')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA core not found"


    # 2. Check for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"
        
    # 3. Analyze fatty acid chain
    # Get the carbon adjacent to the carbonyl in the thioester
    for match in thioester_matches:
        carbonyl_c_index = match[0] # Index of the carbon in C=O
        thio_s_index = match[1]   # Index of the sulfur in -S-

        
        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_index)
        
        #Get the carbon connected to the carbonyl carbon
        connected_carbon_indices = [neighbor.GetIdx() for neighbor in carbonyl_c_atom.GetNeighbors() if neighbor.GetIdx() != thio_s_index and neighbor.GetAtomicNum() == 6]
        if len(connected_carbon_indices) != 1:
             return False, "Could not identify the connecting carbon to the carbonyl in thioester."
        
        start_carbon_index = connected_carbon_indices[0]


        # Count the carbons in the chain
        chain_length = 0
        visited_atoms = set()
        current_atom_index = start_carbon_index
        prev_atom_index = carbonyl_c_index

        while True:
           
            chain_length += 1
            visited_atoms.add(current_atom_index)
            current_atom = mol.GetAtomWithIdx(current_atom_index)
           
           # Find the next carbon in the chain (that is not previous or thioester C=O)
            next_carbon_index = None
            for neighbor in current_atom.GetNeighbors():
                neighbor_index = neighbor.GetIdx()
                if neighbor.GetAtomicNum() == 6 and neighbor_index != prev_atom_index and neighbor_index not in visited_atoms:
                    next_carbon_index = neighbor_index
                    break

            if next_carbon_index is None:
                break #End of chain
            
            prev_atom_index = current_atom_index
            current_atom_index = next_carbon_index

        if  not 6 <= chain_length <= 12:
             return False, f"Fatty acid chain has {chain_length} carbons, should be 6-12."
        
        # Allow for possible double bonds in the chain (done previously by smarts)
        # Check for substituents
        # We do not check here for any substituens.

    return True, "Medium-chain fatty acyl-CoA identified"
"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: CHEBI:26666 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import HydrogenAdder, rdmolops, rdchem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is an unbranched, fully saturated carbon chain ending with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to accurately count atoms
    mol = Chem.AddHs(mol)

    # Check for carboxylic acid group (-C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify the carboxylic acid carbon
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylic_carbons = [match[0] for match in matches]  # First carbon in the pattern
    
    # Generate the molecular graph
    mol_graph = rdmolops.GetAdjacencyMatrix(mol)
    
    # Check for cycles (acyclic)
    if Chem.GetSSSR(mol) > 0:
        return False, "Molecule contains rings"

    # Find all carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found"

    # Check for unsaturation (double/triple bonds)
    unsaturated_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
            # Allow double bond in carboxylic acid group
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            if begin_atom in carboxylic_carbons or end_atom in carboxylic_carbons:
                continue
            unsaturated_bond = True
            break
    if unsaturated_bond:
        return False, "Unsaturated bonds found in carbon chain"

    # Identify the longest carbon chain starting from carboxylic acid carbon
    from collections import deque

    max_chain_length = 0
    is_straight_chain = True

    for carboxylic_carbon in carboxylic_carbons:
        visited = set()
        queue = deque()
        queue.append((carboxylic_carbon, 0, -1))  # (current_atom, chain_length, previous_atom)
        while queue:
            current_atom, chain_length, prev_atom = queue.popleft()
            if current_atom in visited:
                continue
            visited.add(current_atom)
            atom = mol.GetAtomWithIdx(current_atom)
            if atom.GetAtomicNum() != 6 and current_atom != carboxylic_carbon:
                continue  # Skip non-carbon atoms except for the carboxylic carbon
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() != prev_atom]
            if len(neighbors) > 1:
                # More than two neighbors means branching (except for carboxylic carbon)
                if current_atom != carboxylic_carbon:
                    is_straight_chain = False
                    break
            for neighbor in neighbors:
                queue.append((neighbor, chain_length + 1, current_atom))
            if chain_length > max_chain_length:
                max_chain_length = chain_length
        if not is_straight_chain:
            break

    if not is_straight_chain:
        return False, "Carbon chain is branched"

    if max_chain_length < 2:
        return False, "Carbon chain is too short"

    return True, "Molecule is a straight-chain saturated fatty acid"

# Test the function with an example
# smiles_example = 'CCCCCCCCCCCCCCCCCC(O)=O'  # Hexadecanoic acid
# result, reason = is_straight_chain_saturated_fatty_acid(smiles_example)
# print(result, reason)
    

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26666',
                              'name': 'straight-chain saturated fatty acid',
                              'definition': 'Any saturated fatty acid lacking a side-chain.',
                              'parents': ['CHEBI:15841']},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}
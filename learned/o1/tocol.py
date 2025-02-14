"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: CHEBI:26416 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton that is substituted at position 2
    by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define chromanol core SMARTS pattern (chroman-6-ol)
    chromanol_smarts = "Oc1ccc2OCCc2c1"
    chromanol_pattern = Chem.MolFromSmarts(chromanol_smarts)
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chromanol core (chroman-6-ol skeleton) found"
    
    # Find the chromanol core match
    chromanol_match = mol.GetSubstructMatch(chromanol_pattern)
    chromanol_atoms = set(chromanol_match)
    
    # Identify substitution at position 2 (atom index 8 in chromanol_smarts)
    position_2_idx = chromanol_match[8]  # Atom at position 2 in the molecule
    
    # Get substituents at position 2
    position_2_atom = mol.GetAtomWithIdx(position_2_idx)
    substituents = [nbr for nbr in position_2_atom.GetNeighbors() if nbr.GetIdx() not in chromanol_atoms]
    if not substituents:
        return False, "No substitution at position 2 of chromanol core"
    
    # Assume the first substituent is the side chain
    side_chain_atom = substituents[0]
    
    # Use BFS to get all atoms in the side chain
    visited = set()
    to_visit = [side_chain_atom.GetIdx()]
    side_chain_atoms = set()
    while to_visit:
        current_idx = to_visit.pop()
        if current_idx not in visited and current_idx not in chromanol_atoms:
            visited.add(current_idx)
            side_chain_atoms.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            for nbr in current_atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in chromanol_atoms:
                    to_visit.append(nbr.GetIdx())
    
    # Check that side chain consists of only carbon and hydrogen atoms
    side_chain_elements = set(mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in side_chain_atoms)
    if not side_chain_elements.issubset({6, 1}):
        return False, "Side chain contains atoms other than carbon and hydrogen"
    
    # Count number of carbons in side chain
    carbon_count = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if carbon_count != 15:
        return False, f"Side chain has {carbon_count} carbons, expected 15 (three isoprenoid units)"
    
    # Count number of double bonds in side chain
    double_bond_count = 0
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in side_chain_atoms and end_idx in side_chain_atoms:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1
    
    if double_bond_count not in [0, 3]:
        return False, f"Side chain has {double_bond_count} double bonds, expected 0 or 3"
    
    return True, "Molecule is a tocol with chromanol core and appropriate side chain"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26416',
                              'name': 'tocol',
                              'definition': 'A chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.',
                              'parents': []},
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
        'stdout': None,
        'num_true_positives': 150,
        'num_false_positives': 4,
        'num_true_negatives': 182407,
        'num_false_negatives': 23,
        'num_negatives': None,
        'precision': 0.974025974025974,
        'recall': 0.8670520231213873,
        'f1': 0.9174311926605504,
        'accuracy': 0.9998521228585199}
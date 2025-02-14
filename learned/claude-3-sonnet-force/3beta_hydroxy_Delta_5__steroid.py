"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: CHEBI:37991 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is a steroid that has a hydroxyl group at the 3 position,
    and a double bond between carbons 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@@]12[C@H]([C@@]([C@H](C1)C)([C@@]2([H])[H])[H])[H]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for 3-hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    has_3_hydroxyl = any(mol.GetAtomWithIdx(idx).GetIsotope() == 3 for match in hydroxyl_matches for idx in match)
    if not has_3_hydroxyl:
        return False, "No hydroxyl group at position 3"
    
    # Look for double bond between C5 and C6
    c5_c6_double_bond_pattern = Chem.MolFromSmarts("[C@H]=[C@@H]")
    c5_c6_double_bond_matches = mol.GetSubstructMatches(c5_c6_double_bond_pattern)
    has_c5_c6_double_bond = any(
        mol.GetAtomWithIdx(idx1).GetIsotope() == 5 and mol.GetAtomWithIdx(idx2).GetIsotope() == 6
        for match in c5_c6_double_bond_matches
        for idx1, idx2 in [(match[0], match[1]), (match[1], match[0])]
    )
    if not has_c5_c6_double_bond:
        return False, "No double bond between C5 and C6"
    
    return True, "Molecule contains a steroid backbone with a 3-hydroxyl group and a double bond between C5 and C6"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:37991',
        'name': '3beta-hydroxy-Delta(5)-steroid',
        'definition': 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.',
        'parents': ['CHEBI:51411', 'CHEBI:35694']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 100,
    'num_false_positives': 3,
    'num_true_negatives': 182291,
    'num_false_negatives': 15,
    'num_negatives': None,
    'precision': 0.9709302325581395,
    'recall': 0.8695652173913043,
    'f1': 0.9166666666666666,
    'accuracy': 0.9998992490118577
}
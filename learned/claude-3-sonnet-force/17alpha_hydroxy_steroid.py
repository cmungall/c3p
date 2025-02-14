"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:63032 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, FragmentMatcher

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is the alpha-stereoisomer of 17-hydroxy steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid scaffold
    steroid_pattern = Chem.MolFromSmarts("[C@]12CC[C@H]3[C@@H]([C@@H]1[C@H](C2)C)CCC4=CC(=O)CC[C@]34C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid scaffold"

    # Check for 17-hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)[C@]12CC[C@H]3[C@@H]([C@@H]1[C@H](C2)C)CCC4=CC(=O)CC[C@]34C")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Does not have a 17-hydroxyl group"

    # Check for alpha stereochemistry at C17
    matcher = FragmentMatcher.FragmentMatcher()
    matcher.addFragmentSmarts("[C@@H](O)[C@@]12CC[C@H]3[C@@H]([C@@H]1[C@H](C2)C)CCC4=CC(=O)CC[C@]34C", True)
    if not matcher.countMatches(mol):
        return False, "Incorrect stereochemistry at C17"

    # Check for steroid properties
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:
        return False, "Too few rings for a steroid"

    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings > 1:
        return False, "Too many aromatic rings for a steroid"

    return True, "Contains steroid scaffold with alpha-configured 17-hydroxyl group"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:63032',
        'name': '17alpha-hydroxy steroid',
        'definition': 'The alpha-stereoisomer of 17-hydroxy steroid.',
        'parents': ['CHEBI:35699', 'CHEBI:51597']
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
    'num_true_positives': 62,
    'num_false_positives': 0,
    'num_true_negatives': 182422,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}
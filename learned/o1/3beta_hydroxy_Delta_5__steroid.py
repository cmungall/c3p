"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid

Definition: 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.'
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6 (Delta5).

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
    
    # Check for steroid backbone (rings A/B/C/D fused)
    steroid_pattern = Chem.MolFromSmarts('\
    [#6]12CC3C(C1)CCC3C4CCC(C2)C4')  # Simplified steroid nucleus pattern
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 3beta-hydroxy group with correct stereochemistry
    hydroxy_3beta_pattern = Chem.MolFromSmarts('\
    [C@@H]1([O])[C@@H](CC2)CC[C@]2(C)[C@@H]1C')  # Simplified pattern
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "No 3beta-hydroxy group with beta orientation found"
    
    # Check for double bond between positions 5 and 6
    delta5_double_bond_pattern = Chem.MolFromSmarts('\
    C1=CC[C@H](C)CC1')  # Simplified pattern for C5=C6 double bond
    if not mol.HasSubstructMatch(delta5_double_bond_pattern):
        return False, "No Delta(5) double bond found between positions 5 and 6"
    
    return True, "Molecule is a 3beta-hydroxy-Delta(5)-steroid"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '3beta-hydroxy-Delta(5)-steroid',
        'definition': 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.'
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.9740,
    'recall': 0.8671,
    'f1': 0.9174,
    'accuracy': 0.9999
}
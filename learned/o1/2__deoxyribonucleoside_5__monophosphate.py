"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:27150 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside 5'-monophosphate consists of a 2'-deoxyribose sugar attached to a nucleobase
    at the 1' position and a phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # 2'-deoxyribose sugar with attachment points at 1', 3', and 5'
    deoxyribose_pattern = Chem.MolFromSmarts("""
        [C@@H]1([O][C@H](C[C@H]1O)CO)
    """)
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Nucleobase attached at 1' position
    nucleobase_attachment = Chem.MolFromSmarts("""
        [C@@H]1([O][C@H](C[C@H]1O)CO)[N,O,S]*
    """)
    if not mol.HasSubstructMatch(nucleobase_attachment):
        return False, "No nucleobase attached at 1' position"

    # Phosphate group at 5' position
    phosphate_pattern = Chem.MolFromSmarts("""
        [C@H](CO[P](=O)([O-])O)O
    """)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group at 5' position"

    # Check absence of hydroxyl group at 2' position (deoxy sugar)
    hydroxyl_2prime = Chem.MolFromSmarts("""
        [C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO
    """)
    if mol.HasSubstructMatch(hydroxyl_2prime):
        return False, "Hydroxyl group found at 2' position, not a deoxyribose"

    return True, "Molecule is a 2'-deoxyribonucleoside 5'-monophosphate"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27150',
        'name': "2'-deoxyribonucleoside 5'-monophosphate",
        'definition': "A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.",
        'parents': ['CHEBI:36976', 'CHEBI:7228']
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
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}
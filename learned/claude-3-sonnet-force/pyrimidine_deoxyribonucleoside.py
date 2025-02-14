"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: CHEBI:35808 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is a deoxyribonucleoside containing a pyrimidine base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyrimidine ring substructure
    pyrimidine_pattern = Chem.MolFromSmarts("n1cccnc1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine ring found"
    
    # Look for deoxyribose substructure
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]([C@H]([C@@H]1[C@H](CO)O[C@@H]1O)O)O")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"
    
    # Check if pyrimidine and deoxyribose are connected
    atom_sets = [set(mol.GetAtomWithIdx(match_idx).GetNeighbors()) for match_idx in mol.GetSubstructMatches(pyrimidine_pattern)[0]]
    deoxyribose_atoms = set(mol.GetSubstructMatch(deoxyribose_pattern))
    if not any(atom_set & deoxyribose_atoms for atom_set in atom_sets):
        return False, "Pyrimidine and deoxyribose not connected"
    
    return True, "Contains a pyrimidine base connected to a deoxyribose sugar"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35808',
        'name': 'pyrimidine deoxyribonucleoside',
        'definition': 'A deoxyribonucleoside containing a pyrimidine base.'
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
    'num_true_positives': 43,
    'num_false_positives': 0,
    'num_true_negatives': 182443,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 0.9997639944286905
}
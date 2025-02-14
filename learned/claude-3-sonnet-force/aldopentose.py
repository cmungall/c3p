"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:18277 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a (potential) aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly 5 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 5:
        return False, f"Found {num_carbons} carbon atoms, aldopentoses must have 5"
    
    # Check for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    # Check for 4 hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 4:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, aldopentoses must have 4"
    
    # Check for ring structures (furanose or pyranose forms)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No ring found, aldopentoses typically have furanose or pyranose forms"
    
    return True, "Contains 5 carbon atoms, 4 hydroxyl groups, an aldehyde group, and a ring structure"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:18277',
        'name': 'aldopentose',
        'definition': 'A pentose with a (potential) aldehyde group at one end.',
        'parents': ['CHEBI:18021', 'CHEBI:16646']
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
    'num_true_positives': 352,
    'num_false_positives': 4,
    'num_true_negatives': 182411,
    'num_false_negatives': 15,
    'num_negatives': None,
    'precision': 0.9888888888888889,
    'recall': 0.9592233009708737,
    'f1': 0.9737533699847448,
    'accuracy': 0.9998932789627203
}
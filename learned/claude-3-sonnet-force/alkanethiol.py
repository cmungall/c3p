"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:25710 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for -SH group
    sh_pattern = Chem.MolFromSmarts("[SH]")
    if not mol.HasSubstructMatch(sh_pattern):
        return False, "No sulfanyl (-SH) group found"
    
    # Look for alkyl chains attached to sulfur
    alkyl_pattern = Chem.MolFromSmarts("[S]~[CX4]")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No alkyl chain attached to sulfur"

    # Check that all atoms are C, H, or S
    valid_atoms = (6, 1, 16)  # C, H, S
    if any(atom.GetAtomicNum() not in valid_atoms for atom in mol.GetAtoms()):
        return False, "Contains atoms other than C, H, S"

    return True, "Contains sulfanyl (-SH) group attached to an alkyl chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25710',
                          'name': 'alkanethiol',
                          'definition': 'An alkanethiol is a compound in which '
                                        'a sulfanyl group, -SH, is attached '
                                        'to an alkyl group.',
                          'parents': ['CHEBI:33832', 'CHEBI:68680']},
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
    'num_true_positives': 229,
    'num_false_positives': 5,
    'num_true_negatives': 182421,
    'num_false_negatives': 29,
    'num_negatives': None,
    'precision': 0.9786324786324787,
    'recall': 0.8875968992248062,
    'f1': 0.9311505157698182,
    'accuracy': 0.9997838377007552}
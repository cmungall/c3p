"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI:36603 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are sphinganine, its homologs and stereoisomers, and their hydroxy and unsaturated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for sphingoid backbone pattern: C-C-C-N-C-C-O
    backbone_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][NH][CH][CH][OH]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No sphingoid backbone found"
    
    # Check for long alkyl chain (>=14 carbons)
    chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "Alkyl chain too short for sphingoid"
    
    # Check for hydroxyl group(s)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"
    
    # Check for optional unsaturation
    unsaturated_pattern = Chem.MolFromSmarts("[CH2][CH][CH][CH]")
    unsaturated_matches = mol.GetSubstructMatches(unsaturated_pattern)
    
    # Count carbons and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 16 or n_count != 1:
        return False, "Incorrect atom counts for sphingoid"
    
    reason = "Contains sphingoid backbone with long alkyl chain and hydroxyl group"
    if unsaturated_matches:
        reason += ", potentially unsaturated"
    
    return True, reason


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36603',
        'name': 'sphingoid',
        'definition': 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.',
        'parents': ['CHEBI:35621', 'CHEBI:36963']
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
    'num_true_positives': 42,
    'num_false_positives': 1,
    'num_true_negatives': 182418,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9767441860465116,
    'recall': 1.0,
    'f1': 0.988235294117647,
    'accuracy': 0.9999458240233187
}
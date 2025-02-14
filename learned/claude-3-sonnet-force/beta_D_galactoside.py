"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:28577 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is a D-galactoside with beta configuration at its anomeric center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for D-galactose substructure
    galactose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@H]([C@H](O1)O)O)O)O)O")
    galactose_matches = mol.GetSubstructMatches(galactose_pattern)
    if not galactose_matches:
        return False, "No D-galactose substructure found"
    
    # Check for beta configuration at anomeric center
    beta_pattern = Chem.MolFromSmarts("[C@H](O)[C@H]1[C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O")
    beta_matches = mol.GetSubstructMatches(beta_pattern)
    if not beta_matches:
        return False, "Not in beta configuration at anomeric center"
    
    # Check for glycosidic bond (-O- linker)
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic bond found"
    
    return True, "Contains D-galactose in beta configuration at the anomeric center and a glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28577',
                          'name': 'beta-D-galactoside',
                          'definition': 'Any D-galactoside having beta-configuration at its anomeric centre.',
                          'parents': ['CHEBI:36432', 'CHEBI:87711']},
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
    'num_true_positives': 185,
    'num_false_positives': 0,
    'num_true_negatives': 186783,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0}
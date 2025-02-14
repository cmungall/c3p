"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:27815 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, FragmentMatcher

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a glucoside in which the glycoside group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for D-glucose substructure
    d_glucose_pattern = Chem.MolFromSmiles("OC[C@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)CO")
    d_glucose_matcher = FragmentMatcher.FragmentMatcher()
    d_glucose_matcher.AddFragmentSmarts(Chem.MolToSmarts(d_glucose_pattern))
    matches = d_glucose_matcher.GetMatchingFragments(mol)
    if not matches:
        return False, "No D-glucose substructure found"
    
    # Look for glycosidic bonds (O-C-O)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CR][OX2]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bond_matches:
        return False, "No glycosidic bonds found"
    
    # Check molecular weight - glucosides typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for glucoside"

    return True, "Contains D-glucose substructure and glycosidic bonds"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27815',
        'name': 'D-glucoside',
        'definition': 'Any glucoside in which the glycoside group is derived from D-glucose.',
        'parents': ['CHEBI:24005']
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
    'num_true_positives': 512,
    'num_false_positives': 108,
    'num_true_negatives': 182303,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8257575757575758,
    'recall': 1.0,
    'f1': 0.9048557997656454,
    'accuracy': 0.9994097335240181
}
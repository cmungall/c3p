"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:35482 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hopane skeleton pattern
    hopane_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@]2([C@H]1CC[C@]1([C@@H]2CC[C@@H]2[C@@]3(C)CC[C@H](C)[C@](C)(C)[C@@H]3CC[C@@]12C)C)C")
    if not mol.HasSubstructMatch(hopane_pattern):
        return False, "No hopane skeleton found"

    # Check molecular weight - hopanoids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for hopanoid"
    
    # Check for presence of rings
    n_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, "Too few rings for hopanoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:
        return False, "Too few carbons for hopanoid"
    if o_count > 5:
        return False, "Too many oxygens for hopanoid"

    return True, "Contains hopane skeleton with appropriate molecular weight and ring structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35482',
        'name': 'hopanoid',
        'definition': 'A triterpenoid based on a  hopane skeleton.',
        'parents': ['CHEBI:35481', 'CHEBI:38258']
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
    'num_true_positives': 237,
    'num_false_positives': 6,
    'num_true_negatives': 177417,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.9751031210907041,
    'recall': 0.9794891183278957,
    'f1': 0.9772929936305732,
    'accuracy': 0.9998750141729957
}
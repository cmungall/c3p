"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: CHEBI:36721 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid conjugated to a functional group that gives additional hydrophilicity or charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for bile acid backbone pattern
    bile_acid_pattern = Chem.MolFromSmarts("[C@H]1[C@]2([C@]3([C@@]([C@@]4([C@]([C@@H](CC4)C)[H])C)(CC2)[H])[H])(CC1)[H]"
    if not mol.HasSubstructMatch(bile_acid_pattern):
        return False, "No bile acid backbone found"
    
    # Look for conjugated groups
    conjugated_groups = ["NCC(O)=O", # Glycine
                         "NCCS(O)(=O)=O", # Taurine
                         "OS(O)(=O)=O", # Sulfate
                         "OC(=O)C(O)C(O)C(O)O", # Glucuronate
                         "OC1OC(CO)C(O)C(O)C1O", # Glucose
                         "NCCC(O)(O)=O", # Coenzyme A
                         ]
    conjugated_pattern = Chem.MolFromSmarts(f"({')('.join(conjugated_groups)})")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if not conjugated_matches:
        return False, "No conjugated groups found"
    
    # Additional checks
    molwt = rdMolDescriptors.CalcExactMolWt(mol)
    if molwt < 400:
        return False, "Molecular weight too low for bile acid conjugate"
    
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 4:
        return False, "Too few oxygens for bile acid conjugate"
    
    return True, "Contains bile acid backbone with conjugated group"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36721',
        'name': 'bile acid conjugate',
        'definition': 'Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule. Molecules used for conjugation are: glycine, taurine (and other amino acids); sulfuric acid (for which the term ''sulfate'' may be used); glucuronic acid (for which the term ''glucuronate'' may be used); glucose and other uncharged sugars; and coenzyme A.',
        'parents': ['CHEBI:15586', 'CHEBI:76605']
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
    'num_true_positives': 30,
    'num_false_positives': 0,
    'num_true_negatives': 182380,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}
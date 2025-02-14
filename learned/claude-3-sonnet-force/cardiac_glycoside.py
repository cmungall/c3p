"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: CHEBI:28807 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside is a steroid lactone containing sugar residues that acts on the contractile force of cardiac muscles.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@]12CCC[C@]3([C@H](CC[C@]13)O)C4=CC(=O)CC[C@]4([H])C2")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for lactone ring (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C1CC(=O)OC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Look for sugar residues
    sugar_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar residues found"
    
    # Count rotatable bonds to verify presence of sugar chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds for sugar chains"
    
    # Check molecular weight - cardiac glycosides typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for cardiac glycoside"
    
    # Count oxygen atoms - should have several from sugar residues
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 8:
        return False, "Too few oxygen atoms for cardiac glycoside"
    
    return True, "Contains steroid backbone with lactone ring and sugar residues"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28807',
        'name': 'cardiac glycoside',
        'definition': 'Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.',
        'parents': ['CHEBI:35489', 'CHEBI:38621']
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
    'num_true_positives': 194,
    'num_false_positives': 4,
    'num_true_negatives': 182349,
    'num_false_negatives': 37,
    'num_negatives': None,
    'precision': 0.9798477157360406,
    'recall': 0.8397386759581882,
    'f1': 0.9062264150943396,
    'accuracy': 0.9997790160865186
}
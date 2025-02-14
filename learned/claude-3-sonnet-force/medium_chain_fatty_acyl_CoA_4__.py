"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:63128 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA oxoanion resulting from deprotonation 
    of the phosphate and diphosphate groups of any medium-chain fatty acyl-CoA, with a chain
    length typically between 6 and 12 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)[C])=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Look for acyl chain (typically between 6-12 carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3]~[CX3]~[CX3]~[CX3]~[CX3]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "No acyl chain found"
    for match in acyl_chain_matches:
        chain_length = sum(1 for atom in mol.GetAtomWithIdx(idx).GetNeighbors() if atom.GetAtomicNum() == 6 for idx in match)
        if 6 <= chain_length <= 12:
            break
    else:
        return False, "Acyl chain length not in medium-chain range (6-12 carbons)"
    
    # Check for deprotonated phosphate groups
    if sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1 and atom.GetAtomicNum() == 8) != 4:
        return False, "Incorrect number of deprotonated phosphate groups"
    
    return True, "Contains CoA backbone with medium-chain acyl group and deprotonated phosphate groups"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:63128',
        'name': 'medium-chain fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups of any medium-chain fatty acyl-CoA; major species at pH 7.3.',
        'parents': ['CHEBI:63143', 'CHEBI:49206']
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
    'num_true_positives': 404,
    'num_false_positives': 15,
    'num_true_negatives': 182392,
    'num_false_negatives': 9,
    'num_negatives': None,
    'precision': 0.9643284379085504,
    'recall': 0.9781420765027322,
    'f1': 0.971154522292994,
    'accuracy': 0.9992458720942091
}
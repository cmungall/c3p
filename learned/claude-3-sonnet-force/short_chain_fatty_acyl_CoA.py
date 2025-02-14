"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:35681 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA has a fatty acid chain of 6 carbons or fewer attached to a coenzyme A (CoA) moiety via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(C)(CO)C(C(=O)NCCC(=O)NCCS)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A (CoA) moiety found"
    
    # Look for fatty acid chain attached via thioester
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])SCC")
    thioester_match = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_match:
        return False, "No thioester bond found to connect fatty acid chain"
    
    # Count carbons in fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_match = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_match:
        return False, "No fatty acid chain found"
    chain_length = max(len(match) for match in fatty_acid_match)
    if chain_length > 6:
        return False, "Fatty acid chain is too long (>6 carbons)"
    
    # Check for additional functional groups on fatty acid chain
    func_group_pattern = Chem.MolFromSmarts("[NX3,OX2,SX2,PX3]")
    func_group_match = mol.GetSubstructMatches(func_group_pattern)
    if func_group_match:
        return False, "Fatty acid chain contains additional functional groups"
    
    return True, "Contains a fatty acid chain of 6 carbons or fewer attached to coenzyme A via a thioester bond"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35681',
        'name': 'short-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal '
                      'condensation of the thiol group of coenzyme A with the '
                      'carboxy group of any short-chain fatty acid.',
        'parents': ['CHEBI:35679']
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
    'num_true_positives': 91,
    'num_false_positives': 5,
    'num_true_negatives': 182406,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.9481481481481481,
    'recall': 0.9786324786324787,
    'f1': 0.9631908165324309,
    'accuracy': 0.9998576043256688
}
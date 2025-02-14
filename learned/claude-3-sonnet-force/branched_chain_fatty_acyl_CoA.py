"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60444 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA results from the formal condensation of the thiol group of
    coenzyme A with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12)C(=O)NCCCC(=O)NCCSC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found"
    
    # Look for branched fatty acid chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3]([CX4,CX3])([CX4,CX3])[CX3,CX2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No branched fatty acid chain found"
    
    # Check for ester linkage between CoA and fatty acid
    ester_pattern = Chem.MolFromSmarts("CCC(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chain too short to be a fatty acid"
    
    # Check molecular weight - typically >800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for branched-chain fatty acyl-CoA"
    
    return True, "Contains coenzyme A backbone with a branched-chain fatty acid attached via an ester bond"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:60444',
        'name': 'branched-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any branched-chain fatty acid.',
        'parents': ['CHEBI:57363', 'CHEBI:57787']
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
    'num_true_positives': 28,
    'num_false_positives': 1,
    'num_true_negatives': 182408,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9655172413793104,
    'recall': 1.0,
    'f1': 0.9825174825174826,
    'accuracy': 0.9999954584201931
}
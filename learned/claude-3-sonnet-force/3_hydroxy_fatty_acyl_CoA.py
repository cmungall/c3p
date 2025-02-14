"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:26076 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for coenzyme A substructure
    coa_pattern = Chem.MolFromSmarts("OP(OP(OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)(=O)O)(=O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain coenzyme A substructure"
    
    # Look for 3-hydroxy fatty acyl chain
    hydroxy_pattern = Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H2][CX4H2]([CX4H2])[CX4H2][CX4H2][CX4H2][CX4H2]([CX3](=[OX1])[CH3])[OX2H,OX1-]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "Does not contain 3-hydroxy fatty acyl chain"
    
    # Check carbonyl connection between CoA and fatty acid
    coa_atoms = mol.GetSubstructMatches(coa_pattern)[0]
    coa_carbonyl = coa_atoms[0]
    fatty_acid_atoms = hydroxy_matches[0]
    fatty_acid_carbonyl = fatty_acid_atoms[9]
    if mol.GetBondBetweenAtoms(coa_carbonyl, fatty_acid_carbonyl).GetBondType() != Chem.BondType.SINGLE:
        return False, "Coenzyme A and fatty acid not connected via carbonyl"
    
    # Count atoms, look for expected range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 25 or o_count < 12:
        return False, "Atom counts too low for 3-hydroxy fatty acyl-CoA"
    
    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800 or mol_wt > 1200:
        return False, "Molecular weight outside expected range for 3-hydroxy fatty acyl-CoA"
    
    return True, "Contains coenzyme A and 3-hydroxy fatty acyl chain connected via carbonyl"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26076',
        'name': '3-hydroxy fatty acyl-CoA',
        'definition': 'A hydroxy fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.',
        'parents': ['CHEBI:36034', 'CHEBI:35581', 'CHEBI:35549']
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
    'num_true_positives': 160,
    'num_false_positives': 5,
    'num_true_negatives': 182417,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.9697265625,
    'recall': 0.9876543209876543,
    'f1': 0.9785723049707534,
    'accuracy': 0.9998786171784357
}
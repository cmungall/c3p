"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:84947 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    This subclass is defined as any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA(4-) substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA(4-) substructure"
    
    # Look for fatty acyl chain
    fatty_acyl_pattern = Chem.MolFromSmarts("CCC(=O)")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    if not fatty_acyl_matches:
        return False, "No fatty acyl chain found"
    
    # Check for saturated 11-12 bond
    for match in fatty_acyl_matches:
        # Find 11th and 12th carbon atoms
        carbon_11 = mol.GetAtomWithIdx(match[0] + 10)
        carbon_12 = mol.GetAtomWithIdx(match[0] + 11)
        
        # Check bond order between them
        bond = mol.GetBondBetweenAtoms(carbon_11.GetIdx(), carbon_12.GetIdx())
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return False, "Found unsaturated bond between carbons 11 and 12"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Fatty acyl chain too short"
    
    return True, "Contains CoA(4-) substructure with saturated bond between carbons 11 and 12 in fatty acyl chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:84947',
        'name': '11,12-saturated fatty acyl-CoA(4-)',
        'definition': 'Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated. This is needed for a reaction together with CHEBI:84947.',
        'parents': ['CHEBI:36865', 'CHEBI:15998']
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
    'num_true_positives': 149,
    'num_false_positives': 0,
    'num_true_negatives': 182436,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 0.9999998582187252
}
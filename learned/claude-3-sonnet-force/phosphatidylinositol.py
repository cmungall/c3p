"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: CHEBI:17559 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is a glycerophosphoinositol with one phosphatidyl group esterified
    to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for inositol backbone pattern
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H](O)[C@H]([C@H]([C@@H]1O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol backbone found"
    
    # Look for phosphate group attached to inositol
    phos_pattern = Chem.MolFromSmarts("OP(O)(=O)")
    phos_matches = mol.GetSubstructMatches(phos_pattern)
    if not phos_matches:
        return False, "No phosphate group attached to inositol"
    
    # Look for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"
    
    # Check molecular weight - phosphatidylinositols typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylinositol"
    
    return True, "Contains inositol backbone with phosphate group and glycerol backbone with 2 fatty acid chains"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17559',
                          'name': 'phosphatidylinositol',
                          'definition': 'Any glycerophosphoinositol having one '
                                        'phosphatidyl group esterified to one '
                                        'of the hydroxy groups of inositol.',
                          'parents': ['CHEBI:25135']},
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
    'num_true_positives': 160,
    'num_false_positives': 8,
    'num_true_negatives': 182396,
    'num_false_negatives': 20,
    'num_negatives': None,
    'precision': 0.9523809523809523,
    'recall': 0.8888888888888888,
    'f1': 0.9198312236287845,
    'accuracy': 0.9995629937445096}
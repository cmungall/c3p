"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:26676 lysophosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate, with one fatty acid chain attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X3][CHX3][CH2X3]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for phosphate group (-O-P(=O)(O)O)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, must have exactly 1"

    # Check for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 1:
        return False, f"Must have exactly one fatty acid chain, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Fatty acid chain too short"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12:
        return False, "Too few carbons for lysophosphatidic acid"
    if o_count != 7:
        return False, "Must have exactly 7 oxygens (1 ester, 1 phosphate)"

    return True, "Contains glycerol backbone with 1 fatty acid chain and 1 phosphate group"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26676',
                                         'name': 'lysophosphatidic acid',
                                         'definition': 'Any monoacylglycerol '
                                                       'phosphate obtained by '
                                                       'hydrolytic removal of '
                                                       'one of the two acyl '
                                                       'groups of any '
                                                       'phosphatidic acid or '
                                                       'derivatives therein.',
                                         'parents': ['CHEBI:36081', 'CHEBI:48246']},
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
                   'num_true_positives': 44,
                   'num_false_positives': 13,
                   'num_true_negatives': 181579,
                   'num_false_negatives': 0,
                   'num_negatives': None,
                   'precision': 0.7718463301205284,
                   'recall': 1.0,
                   'f1': 0.8721804511278195,
                   'accuracy': 0.9997124539221404}
"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:18299 1-phosphatidyl-1D-myo-inositol

A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer and the phosphatidyl group is located at its position 1.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol backbone pattern
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol backbone found"

    # Look for phosphate group attached to position 1
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])[O;X2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 1"

    phosphate_atom = mol.GetAtomWithIdx(list(phosphate_matches[0])[1])
    if phosphate_atom.GetTotalNumHs() != 0:
        return False, "Phosphate group not attached to inositol backbone"

    # Look for two fatty acid chains (long carbon chains attached to glycerol backbone)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) != 1:
        return False, f"Found {len(glycerol_matches)} glycerol backbones, need exactly 1"

    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    return True, "Contains inositol backbone with phosphate group at position 1 and two fatty acid chains"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18299',
                          'name': '1-phosphatidyl-1D-myo-inositol',
                          'definition': 'A phosphatidylinositol in which the '
                                        'inositol moiety is the 1D-myo isomer '
                                        'and the phosphatidyl group is located '
                                        'at its position 1.',
                          'parents': ['CHEBI:17834', 'CHEBI:18291']},
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
    'num_true_positives': 768,
    'num_false_positives': 15,
    'num_true_negatives': 182428,
    'num_false_negatives': 361,
    'num_negatives': None,
    'precision': 0.9808264462809917,
    'recall': 0.6804123711340206,
    'f1': 0.8053097345132744,
    'accuracy': 0.9984106908774476}
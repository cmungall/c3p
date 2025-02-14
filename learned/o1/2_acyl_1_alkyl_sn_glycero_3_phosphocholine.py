"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:63866 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    This class has:
    - A glycerol backbone with sn (stereospecific numbering) configuration at position 2.
    - An ether-linked alkyl chain at position 1.
    - An ester-linked acyl chain at position 2.
    - A phosphocholine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has explicit hydrogens to preserve stereochemistry
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for the required substructures

    # Glycerol backbone with sn stereochemistry at position 2
    # The [C@H] denotes chiral center with specified configuration
    glycerol_sn_pattern = Chem.MolFromSmarts("[CH2][C@H](O*)[CH2]")

    # Ether-linked alkyl chain at position 1 (position attached via oxygen to carbon chain without carbonyl)
    ether_alkyl_pattern = Chem.MolFromSmarts("[CH2]-[O]-[C;H2,H3][C;H2,H3]")

    # Ester-linked acyl chain at position 2 (position attached via oxygen to carbonyl carbon)
    ester_acyl_pattern = Chem.MolFromSmarts("[C@H](O-C(=O)-[C;H2,H3])[CH2]")

    # Phosphocholine group at position 3
    phosphocholine_pattern = Chem.MolFromSmarts("[CH2]-O-P(=O)([O-])-OCC[N+](C)(C)C")

    # Check for glycerol backbone with sn stereochemistry
    if not mol.HasSubstructMatch(glycerol_sn_pattern):
        return False, "No glycerol backbone with sn stereochemistry found"

    # Check for ether-linked alkyl chain at position 1
    if not mol.HasSubstructMatch(ether_alkyl_pattern):
        return False, "No ether-linked alkyl chain at position 1 found"

    # Check for ester-linked acyl chain at position 2
    if not mol.HasSubstructMatch(ester_acyl_pattern):
        return False, "No ester-linked acyl chain at position 2 found"

    # Check for phosphocholine group at position 3
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group at position 3 found"

    # If all patterns match, molecule is classified as 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
    return True, "Molecule matches all required substructures for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:63866',
        'name': '2-acyl-1-alkyl-sn-glycero-3-phosphocholine',
        'definition': 'An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl groups are located at positions 1 and 2 respectively.',
        'parents': ['CHEBI:63865', 'CHEBI:63864']
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}
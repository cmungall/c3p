"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: CHEBI:17680 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a ceramide with a glucose moiety attached to the ceramide backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ceramide backbone pattern
    ceramide_pattern = Chem.MolFromSmarts("[NX3H2;$(NC(=O)C)][CX4H](C[CX4])([CX4])[CX3](=O)[OX2H]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Look for glucose pattern
    glucose_pattern = Chem.MolFromSmarts("[OX2;$([C@H]1[C@@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O[CX4;!$(NC=O)])][CX4]")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No glucose moiety found"

    # Check for ether/ester bond between glucose and ceramide backbone
    ether_pattern = Chem.MolFromSmarts("[OX2;$([C@H]1[C@@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O[CX4;!$(NC=O)])][CX4;$(NC(=O)C)]")
    ester_pattern = Chem.MolFromSmarts("[OX2;$([C@H]1[C@@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O[CX3]=O)]")
    if not (mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No ether/ester bond between glucose and ceramide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for glucosylceramide"
    if o_count < 8:
        return False, "Too few oxygens for glucosylceramide"

    return True, "Contains a ceramide backbone with a glucose moiety attached"

__metadata__ = {
    'chemical_class': {'id': 'CHEBI:17680',
                       'name': 'glucosylceramide',
                       'definition': 'Any of the cerebrosides in which the monosaccharide head group is glucose.',
                       'parents': ['CHEBI:35783', 'CHEBI:83427']},
    'config': {'llm_model_name': 'lbl/claude-sonnet',
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
    'num_true_positives': 24,
    'num_false_positives': 1,
    'num_true_negatives': 182406,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.96,
    'recall': 1.0,
    'f1': 0.9795918367346939,
    'accuracy': 0.9999944795918365}
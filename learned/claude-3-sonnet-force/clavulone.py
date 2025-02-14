"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: CHEBI:137276 clavulone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are esterified prostanoids obtained from marine corals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for prostanoid backbone pattern
    prostanoid_pattern = Chem.MolFromSmarts("[C@@H](CC=CCCC)(CC=C)"
    if not mol.HasSubstructMatch(prostanoid_pattern):
        return False, "No prostanoid backbone found"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Look for halogen substituents (Cl, Br, I)
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    halogen_matches = mol.GetSubstructMatches(halogen_pattern)
    if not halogen_matches:
        return False, "No halogen substituents found"

    # Check for long carbon chains (>5 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]"
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if not long_chain_matches:
        return False, "No long carbon chains found"

    return True, "Contains prostanoid backbone, ester groups, halogen substituents, and long carbon chains"

__metadata__ = {'chemical_class': {'id': 'CHEBI:137276',
                                   'name': 'clavulone',
                                   'definition': 'A class of esterified prostanoids obtained from marine corals.',
                                   'parents': ['CHEBI:23058', 'CHEBI:38513', 'CHEBI:50114']},
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
                 'num_true_positives': 166,
                 'num_false_positives': 14,
                 'num_true_negatives': 182386,
                 'num_false_negatives': 18,
                 'num_negatives': None,
                 'precision': 0.9221556886227545,
                 'recall': 0.9021505376344086,
                 'f1': 0.9119743521594684,
                 'accuracy': 0.9999203358688792}
"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem is a beta-lactam antibiotic with a carbapenem skeleton,
    which is a fused beta-lactam ring and a five-membered ring with an exocyclic double bond,
    variously substituted at positions 3, 4, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carbapenem core structure using SMILES
    # This represents the fused beta-lactam ring and five-membered ring with an exocyclic double bond
    core_smiles = 'C1C=CN2C1C(=O)C2'  # Simplified carbapenem core structure
    core_mol = Chem.MolFromSmiles(core_smiles)
    if core_mol is None:
        return False, "Error in defining carbapenem core structure"

    # Check if molecule contains the carbapenem core
    if not mol.HasSubstructMatch(core_mol):
        return False, "Does not contain carbapenem core structure"

    # If needed, further checks for substitutions at positions 3, 4, and 6 can be added here

    return True, "Contains carbapenem core structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'carbapenems',
        'definition': 'The class of beta-lactam antibiotics whose members have a carbapenem skeleton which is variously substituted at positions 3, 4, and 6.',
        'parents': ['CHEBI:19258', 'CHEBI:26836']
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
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}
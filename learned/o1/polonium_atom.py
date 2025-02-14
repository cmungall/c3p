"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.
    A polonium atom is a single atom with atomic number 84.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polonium atom, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that there is exactly one atom
    num_atoms = mol.GetNumAtoms()
    if num_atoms != 1:
        return False, f"Molecule has {num_atoms} atoms; expected a single atom"

    # Get the atom
    atom = mol.GetAtomWithIdx(0)

    # Check atomic number
    atomic_num = atom.GetAtomicNum()
    if atomic_num != 84:
        return False, f"Atomic number is {atomic_num}; expected 84 for polonium"

    # Check for isotope information
    isotope = atom.GetIsotope()
    if isotope == 0:
        return True, "Polonium atom without specific isotope"

    # List of known polonium isotopes (from examples)
    known_isotopes = {190, 191, 193, 194, 195, 196, 197, 198, 199,
                      200, 201, 202, 203, 204, 206, 208, 209, 210,
                      211, 214, 216, 217, 218}

    if isotope in known_isotopes:
        return True, f"Polonium atom with isotope {isotope}"
    else:
        return False, f"Isotope {isotope} is not a known polonium isotope"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'polonium atom',
        'definition': 'A radioactive metallic element discovered in 1898 by Marie Sklodowska Curie and named after her home country, Poland (Latin Polonia).',
        'parents': []
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
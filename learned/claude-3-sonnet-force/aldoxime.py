"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:35781 aldoxime

An aldoxime is an oxime of an aldehyde, with the general structure RCH=NOH.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is an oxime of an aldehyde, with the general structure RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for oxime functional group (-CH=N-O-)
    oxime_pattern = Chem.MolFromSmarts("[CH]=N[OH]")
    oxime_matches = mol.GetSubstructMatches(oxime_pattern)
    if not oxime_matches:
        return False, "No oxime functional group (-CH=N-O-) found"

    # Check if oxime is attached to an aldehyde (not a ketone)
    for match in oxime_matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetDegree() > 1:
            return True, "Contains aldoxime functional group (-CH=N-O-)"

    return False, "Oxime not attached to aldehyde carbon"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35781',
        'name': 'aldoxime', 
        'definition': 'Oximes of aldehydes RCH=NOH.'
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
    'num_true_positives': 1196,
    'num_false_positives': 15, 
    'num_true_negatives': 182400,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9875515695067265,
    'recall': 1.0,
    'f1': 0.9937106918238994,
    'accuracy': 0.9991828202390213
}
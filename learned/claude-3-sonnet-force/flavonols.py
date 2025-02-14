"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:17773 flavonol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonol(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is a flavonoid (has benzopyran backbone)
    flavonoid_pattern = Chem.MolFromSmarts("[c1]2[c]([c]3[c]([c]([c]2[o]1)[O])[O])-[c]([c]([c]3=O)[O])")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Not a flavonoid (missing benzopyran backbone)"

    # Check for hydroxy group at position 3
    position_3_oh = Chem.MolFromSmarts("[c1]2[c]([c]([c]([o]2)[OH])[O])-[c]([c]([c]1=O)[O])")
    if not mol.HasSubstructMatch(position_3_oh):
        return False, "No hydroxy group at position 3 of heterocyclic ring"

    return True, "Contains a hydroxyflavone backbone with a hydroxy group at position 3"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17773',
        'name': 'flavonol',
        'definition': 'Any hydroxyflavone in which is the ring hydrogen at position 3 of the heterocyclic ring is replaced by a hydroxy group.',
        'parents': ['CHEBI:17282']
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
    'num_true_positives': 176,
    'num_false_positives': 1,
    'num_true_negatives': 182409,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9943502824858757,
    'recall': 1.0,
    'f1': 0.9971547266730862,
    'accuracy': 0.9999455453326621
}
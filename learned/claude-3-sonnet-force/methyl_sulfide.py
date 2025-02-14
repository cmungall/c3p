"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:33315 methyl sulfide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any aliphatic sulfide where at least one of the organyl groups
    attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfur atom
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atom found"
    
    # Check if any sulfur atom has a methyl group attached
    for sulfur in sulfur_atoms:
        neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in sulfur.GetNeighbors()]
        methyl_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 6 and len(neighbor.GetNeighbors()) == 4]
        if methyl_neighbors:
            return True, "Contains a sulfur atom with a methyl group attached"
    
    return False, "No sulfur atom found with a methyl group attached"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33315',
        'name': 'methyl sulfide',
        'definition': 'Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.',
        'parents': ['CHEBI:33308', 'CHEBI:33311']
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
    'num_true_positives': 41,
    'num_false_positives': 0,
    'num_true_negatives': 182470,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 0.9999778087369353
}
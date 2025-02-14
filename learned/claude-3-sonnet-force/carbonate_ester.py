"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: CHEBI:35462 carbonate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester is a carbonic acid derivative where the hydrogens are replaced by organyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbonate ester functional group (-O-C(=O)-O-)
    carbonate_ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[OX2]")
    carbonate_ester_matches = mol.GetSubstructMatches(carbonate_ester_pattern)
    if not carbonate_ester_matches:
        return False, "No carbonate ester functional group found"
    
    # Check for organyl groups (alkyl/aryl groups) attached to carbonate ester
    for match in carbonate_ester_matches:
        carbonate_ester_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        organyl_attached = False
        for atom in carbonate_ester_atoms:
            if atom.GetAtomicNum() == 8:  # Oxygen atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetHybridization() in (Chem.HybridizationType.SP3, Chem.HybridizationType.SP2):
                        organyl_attached = True
                        break
        if organyl_attached:
            return True, "Contains carbonate ester functional group with organyl groups attached"
    
    return False, "No organyl groups attached to carbonate ester"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35462',
        'name': 'carbonate ester',
        'definition': 'Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups.',
        'parents': ['CHEBI:24863', 'CHEBI:50114']
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
    'num_true_positives': 673,
    'num_false_positives': 3,
    'num_true_negatives': 185902,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.9957098957098957,
    'recall': 0.9912528480704371,
    'f1': 0.9934783367951536,
    'accuracy': 0.9999267758896635
}
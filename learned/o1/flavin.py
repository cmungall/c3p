"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: flavin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of the dimethylisoalloxazine skeleton with a substituent on the 10 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the dimethylisoalloxazine core SMARTS pattern
    # This pattern matches the tricyclic ring system with methyl groups at positions 7 and 8
    flavin_core_smarts = 'Cc1cc2nc3c(nc(=O)[nH]c3=O)n(c2cc1C)[#7]'
    flavin_core = Chem.MolFromSmarts(flavin_core_smarts)
    if flavin_core is None:
        return None, "Error in defining flavin core SMARTS pattern"

    # Check for flavin core
    if not mol.HasSubstructMatch(flavin_core):
        return False, "No dimethylisoalloxazine core found"

    # Get the matching substructure
    matches = mol.GetSubstructMatches(flavin_core)
    if not matches:
        return False, "No match for the flavin core"

    # Check for substituent at position 10
    # Position 10 is the nitrogen atom in the isoalloxazine ring system
    # We need to check if this nitrogen has a substituent (i.e., connected to any atom outside the core)
    substituent_found = False
    for match in matches:
        flavin_atoms = set(match)
        # Index of the N atom at position 10 in the SMARTS pattern
        idx_N10 = match[13]
        N10_atom = mol.GetAtomWithIdx(idx_N10)
        for neighbor in N10_atom.GetNeighbors():
            if neighbor.GetIdx() not in flavin_atoms:
                substituent_found = True
                break
        if substituent_found:
            break

    if not substituent_found:
        return False, "No substituent found at position 10"

    return True, "Contains dimethylisoalloxazine core with substituent at position 10"


__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'flavin',
                              'definition': 'A derivative of the dimethylisoalloxazine '
                                            '(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) '
                                            'skeleton, with a substituent on the 10 position.',
                              'parents': []},
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
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}
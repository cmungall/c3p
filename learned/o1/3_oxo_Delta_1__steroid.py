"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:35164 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is any 3-oxo steroid that contains a double bond between positions 1 and 2.
    Due to limitations in assigning atom positions, this function checks for the steroid backbone,
    the presence of at least one ketone group in a ring, and at least one double bond in a ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern (four fused rings)
    steroid_pattern = Chem.MolFromSmarts('[*]12[*]3[*]4[*]1[*]3[*]2[*]4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for ketone groups (C=O)
    ketone_pattern = Chem.MolFromSmarts('C(=O)')
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_in_ring = False
    for match in ketone_matches:
        carbon_idx = match[0]  # Carbon atom in C=O
        if mol.GetAtomWithIdx(carbon_idx).IsInRing():
            ketone_in_ring = True
            break
    if not ketone_in_ring:
        return False, "No ketone group found in ring"

    # Check for double bonds in rings
    double_bond_in_ring = False
    for bond in mol.GetBonds():
        if bond.IsInRing() and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_in_ring = True
            break
    if not double_bond_in_ring:
        return False, "No double bond found in ring"

    return True, "Contains steroid backbone with ketone and double bond in rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35164',
                              'name': '3-oxo-Delta(1) steroid',
                              'definition': 'Any 3-oxo steroid that contains a '
                                            'double bond between positions 1 and 2.'},
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
        'stdout': None}
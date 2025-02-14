"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is a sterol ester obtained by formal condensation of the carboxy group
    of any carboxylic acid with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define cholesterol backbone SMARTS pattern (sterane skeleton with stereochemistry)
    cholesterol_pattern = Chem.MolFromSmarts(
        '[C@@H]1(CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@@]4(C)[C@@H]3CC=C2C1)O)[C@H](C)CCCC(C)C'
    )
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "Cholesterol backbone not found"

    # Define ester bond at 3-position (oxygen connected to C3 and acyl group)
    ester_bond_pattern = Chem.MolFromSmarts(
        '[C@@H]1(CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@@]4(C)[C@@H]3CC=C2C1)O[C]=O)[C@H](C)CCCC(C)C'
    )
    ester_bond_query = Chem.MolFromSmarts('OC(=O)')
    ester_bond_matches = mol.GetSubstructMatches(ester_bond_query)
    if not ester_bond_matches:
        return False, "No ester bond found at position 3"

    # Check that the ester is at position 3
    # Get the atom index of the oxygen in the ester bond
    ester_oxygen = None
    for match in ester_bond_matches:
        o_idx = match[0]  # oxygen atom
        # Check if this oxygen is connected to the cholesterol backbone at C3
        atom = mol.GetAtomWithIdx(o_idx)
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.HasSubstructMatch(cholesterol_pattern):
                ester_oxygen = o_idx
                break
        if ester_oxygen is not None:
            break
    if ester_oxygen is None:
        return False, "Ester bond not connected to cholesterol backbone at position 3"

    # Confirm that the acyl group is derived from a carboxylic acid (R-C(=O)-O-)
    # This can be inferred from the ester pattern already matched

    return True, "Molecule is a cholesteryl ester (cholesterol esterified at position 3)"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cholesteryl ester',
        'definition': 'A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of cholesterol.',
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
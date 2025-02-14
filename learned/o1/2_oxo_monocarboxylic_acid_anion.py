"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:25450 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion is an oxo monocarboxylic acid anion in which the oxo group is located at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylate group (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate = len(carboxylate_matches)
    if num_carboxylate != 1:
        return False, f"Found {num_carboxylate} carboxylate groups, require exactly 1"

    # Get the index of the carboxylate carbon atom
    carboxylate_C_idx = carboxylate_matches[0][0]
    carboxylate_C = mol.GetAtomWithIdx(carboxylate_C_idx)

    # Find the alpha carbon (adjacent to the carboxylate carbon)
    alpha_carbons = [nbr for nbr in carboxylate_C.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(alpha_carbons) == 0:
        return False, "No alpha carbon adjacent to carboxylate group"

    # Check if any alpha carbon has an oxo group (=O) attached
    has_oxo_at_alpha = False
    for alpha_C in alpha_carbons:
        # Look for a double bond to oxygen on the alpha carbon
        oxo_found = False
        for bond in alpha_C.GetBonds():
            neighbor = bond.GetOtherAtom(alpha_C)
            if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                oxo_found = True
                break
        if oxo_found:
            has_oxo_at_alpha = True
            break

    if not has_oxo_at_alpha:
        return False, "No oxo group (=O) at the alpha carbon (2-position)"

    return True, "Molecule is a 2-oxo monocarboxylic acid anion"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25450',
        'name': '2-oxo monocarboxylic acid anion',
        'definition': 'An oxo monocarboxylic acid anion in which the oxo group is located at the 2-position.',
        'parents': ['CHEBI:35756', 'CHEBI:37549']
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
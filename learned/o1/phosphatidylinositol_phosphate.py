"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate is a glycerophospholipid with a glycerol backbone,
    two fatty acid chains, a phosphate group attached to the sn-3 position, and an inositol ring
    attached via the phosphate group. The inositol ring may have additional phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with two esterified fatty acid chains
    # Glycerol backbone pattern with esterified fatty acids at sn-1 and sn-2
    glycerol_fatty_acid_pattern = Chem.MolFromSmarts("""
        [C@@H](CO[C@@H](C(=O)[#6])[#6])[CH2]O[P](=O)(O*)O[C@H]1[C@@H](O*)[C@H](O*)[C@@H](O*)[C@H](O*)[C@H]1O*
    """)
    if not mol.HasSubstructMatch(glycerol_fatty_acid_pattern):
        return False, "No glycerol backbone with two esterified fatty acid chains found"

    # Check for phosphate group attached to sn-3 position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to glycerol backbone"

    # Check for inositol ring attached via phosphate group
    # Inositol ring pattern (6-membered ring with five hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("""
        O[P](=O)(O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O
    """)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring attached via phosphate group"

    # Check for additional phosphate groups on inositol ring
    # Phosphorylated inositol positions (variable number of phosphates)
    phosphate_on_inositol_pattern = Chem.MolFromSmarts("""
        [C@H]1([C@H]([C@@H]([C@H]([C@H]([C@@H]1O[P](=O)(O)O)O[P](=O)(O)O)O[P](=O)(O)O)O)O)O
    """)
    phosphate_count = len(mol.GetSubstructMatches(phosphate_on_inositol_pattern))

    # Count total phosphate groups attached to inositol ring
    inositol_phosphate_positions = [3, 4, 5]
    total_phosphates = 0
    for position in inositol_phosphate_positions:
        phosphate_smarts = f"""
            [C@H]1[C@H]({{O}})[C@@H]({{O}})[C@H]({{O}})[C@H]({{O}})[C@@H]1O[P](=O)(O)O
        """.format(O='O' * position)
        phosphate_position_pattern = Chem.MolFromSmarts(phosphate_smarts)
        if mol.HasSubstructMatch(phosphate_position_pattern):
            total_phosphates += 1

    # Phosphatidylinositol phosphates have up to 3 phosphates on the inositol ring
    if total_phosphates == 0:
        return False, "No additional phosphate groups on inositol ring found"

    # Overall validation passed
    return True, "Molecule is a phosphatidylinositol phosphate"

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'phosphatidylinositol phosphate',
        'definition': 'Any member of the phosphoinositide family of compounds, of which seven occur naturally.',
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
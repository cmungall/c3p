"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:28874 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    This involves checking for a glycerol backbone esterified at positions 1 and 2 with fatty acids,
    a phosphate group at position 3, and a 1D-myo-inositol ring attached to the phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone with ester groups at positions 1 and 2
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO[P](=O)(O)O)(COC(=O)*)OC(=O)*")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct substitution found"
    
    # Check for phosphate group at the third carbon of glycerol
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group linked to glycerol found"
    
    # Check for inositol ring attached to phosphate group
    inositol_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol ring attached to phosphate group found"
    
    # Check for correct stereochemistry in inositol ring (1D-myo-inositol)
    stereo_inositol = Chem.MolFromSmiles("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(stereo_inositol):
        return False, "Inositol ring does not have correct stereochemistry (1D-myo-inositol)"
    
    # Ensure fatty acid chains are attached via ester bonds at positions 1 and 2
    ester_pattern = Chem.MolFromSmarts("COC(=O)[C][C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Less than two esterified fatty acid chains found"
    
    # Optionally, check for long fatty acid chains (common in phosphatidylinositols)
    num_carbon_chains = 0
    for ester in ester_matches:
        # Extract the ester group and count carbons in the fatty acid chain
        ester_atom_idx = ester[0]
        atom = mol.GetAtomWithIdx(ester_atom_idx)
        chain_length = 0
        visited_atoms = set()
        stack = [atom.GetIdx()]
        while stack:
            idx = stack.pop()
            if idx in visited_atoms:
                continue
            visited_atoms.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                chain_length += 1
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if neighbor.GetAtomicNum() == 6 and nbr_idx not in visited_atoms:
                    stack.append(nbr_idx)
        if chain_length >= 8:
            num_carbon_chains += 1
    if num_carbon_chains < 2:
        return False, "Fatty acid chains are too short"
    
    return True, "Molecule is a 1-phosphatidyl-1D-myo-inositol with correct structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28874',
        'name': '1-phosphatidyl-1D-myo-inositol',
        'definition': 'A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer and the phosphatidyl group is located at its position 1.',
        'parents': ['CHEBI:49117', 'CHEBI:24851']
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}
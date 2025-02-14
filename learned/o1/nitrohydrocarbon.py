"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:XXXX nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens has been replaced by nitro groups (C–NO2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that all atoms are C, H, N, or O
    allowed_atomic_nums = {1, 6, 7, 8}  # H, C, N, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom other than C, H, N, and O: {atom.GetSymbol()}"

    # Check for nitro group(s) (–NO2)
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')  # Nitro group pattern
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"

    # Check that nitro groups are attached to carbons
    nitro_carbons = False
    for match in nitro_matches:
        nitro_n_idx = match[0]  # Index of nitrogen in nitro group
        nitro_n = mol.GetAtomWithIdx(nitro_n_idx)
        attached_atoms = nitro_n.GetNeighbors()
        carbon_found = False
        for neighbor in attached_atoms:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                carbon_found = True
            elif neighbor.GetAtomicNum() != 8:  # Should only be attached to O and C
                return False, "Nitro group is not attached to carbon"
        if not carbon_found:
            return False, "Nitro group is not attached to carbon"
        else:
            nitro_carbons = True
    if not nitro_carbons:
        return False, "Nitro groups are not attached to carbon"

    # Remove nitro groups to check the remaining structure
    mol_no_nitro = Chem.RWMol(mol)
    atoms_to_remove = set()
    for match in nitro_matches:
        for idx in match:
            atoms_to_remove.add(idx)
    # Remove atoms in reverse order to keep indices valid
    for idx in sorted(atoms_to_remove, reverse=True):
        mol_no_nitro.RemoveAtom(idx)

    # Check that remaining atoms are only C and H
    for atom in mol_no_nitro.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, f"Contains atom other than C and H in the hydrocarbon part: {atom.GetSymbol()}"

    return True, "Contains hydrocarbon with nitro group(s) attached to carbon"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXX',
        'name': 'nitrohydrocarbon',
        'definition': 'A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.',
        'parents': ['CHEBI:35474', 'CHEBI:62834']
    },
    'config': {
        'llm_model_name': 'INSERT_MODEL_NAME_HERE',
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
"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36737 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as having a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a fully conjugated cyclic dione pattern
    quinone_pattern = Chem.MolFromSmarts("[$(C=C)$(C(=O)C=C),$(C=C)$(C=C(=O))&$(RingMembership)]")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No fully conjugated cyclic dione structure found"

    # Check if the dione is part of a larger conjugated system
    conjugated_system = AllChem.GetConnectedBnds(mol, asMerge=True)
    if len(conjugated_system) < 6:
        return False, "Dione not part of a larger conjugated system"
    
    # Check if the conjugated system is cyclic
    if not any(ring.IsAromaticRing() for ring in AllChem.GetSSSR(mol)):
        return False, "Conjugated system is not cyclic"

    # Check for other common quinone features (optional)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 2:
        return False, "Fewer than 2 oxygens in the molecule"

    n_carbonyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[C=O]")))
    if n_carbonyls < 2:
        return False, "Fewer than 2 carbonyl groups in the molecule"

    return True, "Contains a fully conjugated cyclic dione structure (quinone)"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36737',
        'name': 'quinone',
        'definition': 'Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones, derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups with any necessary rearrangement of double bonds (polycyclic and heterocyclic analogues are included).',
        'parents': ['CHEBI:24631', 'CHEBI:33567']
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
    'num_true_positives': 81,
    'num_false_positives': 3,
    'num_true_negatives': 182688,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.9641691611886184,
    'recall': 0.9310344827586207,
    'f1': 0.9473684210526315,
    'accuracy': 0.9998547739944285
}
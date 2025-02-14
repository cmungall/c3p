"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: CHEBI:47802 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside with a triterpenoid terpene moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for triterpenoid backbone
    triterpenoid_pattern = Chem.MolFromSmarts("[C@@H]1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@]5([C@@](CC4)CCC5)C)[H])=CC3)C)(CC2)C)[H]CCC1")
    if not mol.HasSubstructMatch(triterpenoid_pattern):
        return False, "No triterpenoid backbone found"

    # Look for sugar moieties (any ring with -O- atoms attached)
    sugar_pattern = Chem.MolFromSmarts("[OR1]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moieties found"

    # Count rotatable bonds to verify presence of sugar chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Insufficient rotatable bonds for sugar chains"

    # Check molecular weight - triterpenoid saponins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for triterpenoid saponin"

    # Check elemental composition - triterpenoid saponins contain C, H, O
    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    if set(atoms) - set([6, 1, 8]):
        return False, "Contains elements other than C, H, O"

    return True, "Contains triterpenoid backbone and sugar moieties"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:47802',
        'name': 'triterpenoid saponin',
        'definition': 'A terpene glycoside in which the terpene moiety is a triterpenoid.',
        'parents': ['CHEBI:36334', 'CHEBI:24702']
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
    'num_true_positives': 184,
    'num_false_positives': 24,
    'num_true_negatives': 182382,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.8846153846153846,
    'recall': 0.9892473118279569,
    'f1': 0.9337349397590361,
    'accuracy': 0.9998615447155333
}
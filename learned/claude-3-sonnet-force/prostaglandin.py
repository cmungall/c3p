"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:17915 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring compounds derived from the parent C20 acid, prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for prostanoic acid core structure
    prostanoic_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@H]([C@@]1(CCC(=O)O)[H])[H])[C@@H](C[C@@]2([H])O)C=O")
    if not mol.HasSubstructMatch(prostanoic_pattern):
        return False, "Does not contain prostanoic acid core structure"
    
    # Count carbon atoms - prostaglandins have 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Incorrect number of carbons ({c_count}), prostaglandins have 20"
    
    # Count oxygen atoms - prostaglandins have 3-6 oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3 or o_count > 6:
        return False, f"Incorrect number of oxygens ({o_count}), prostaglandins have 3-6"
    
    # Count rings - prostaglandins have 2 rings
    ring_info = mol.GetRingInfo()
    n_rings = len(set(x for y in ring_info.AtomRings() for x in y))
    if n_rings != 2:
        return False, f"Incorrect number of rings ({n_rings}), prostaglandins have 2"
    
    # Check for trans-double bond in side chain
    side_chain_pattern = Chem.MolFromSmarts("[CX4]/C=C/[CX4]")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "Missing trans-double bond in side chain"
    
    # Check for cis-double bonds in the cyclopentenone ring
    cis_ring_pattern = Chem.MolFromSmarts("[CX3]/C=C/[CX3]")
    cis_ring_matches = mol.GetSubstructMatches(cis_ring_pattern)
    if len(cis_ring_matches) < 1:
        return False, "Missing cis-double bond in cyclopentenone ring"
    
    # All checks pass, classify as prostaglandin
    return True, "Contains prostanoic acid core structure with characteristic ring systems and side chain"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17915',
        'name': 'prostaglandin',
        'definition': 'Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.',
        'parents': ['CHEBI:36594', 'CHEBI:51305']
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
    'num_true_positives': 235,
    'num_false_positives': 3,
    'num_true_negatives': 182374,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9874212535345849,
    'recall': 1.0,
    'f1': 0.9936576614967748,
    'accuracy': 0.9998361457454583
}
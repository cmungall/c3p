"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: CHEBI:35838 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is a saturated fatty acid lacking a side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for straight carbon chain (no branches)
    straight_chain_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if mol.HasSubstructMatch(straight_chain_pattern):
        return False, "Branched carbon chain found"
    
    # Check for saturation (no double or triple bonds)
    if not AllChem.EmbeddedMolecule.GetAtomQuerySet(mol, "aliphaticOnly"):
        return False, "Unsaturated bonds found"
    
    # Check for hydroxy groups (optional, but no other substituents)
    hydroxy_pattern = Chem.MolFromSmarts("O[H]")
    other_substituents = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):  # H, C, O
            other_substituents.append(atom.GetSymbol())
    if other_substituents:
        return False, f"Found other substituents: {', '.join(other_substituents)}"
    
    return True, "Straight-chain saturated fatty acid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35838',
        'name': 'straight-chain saturated fatty acid',
        'definition': 'Any saturated fatty acid lacking a side-chain.',
        'parents': ['CHEBI:35839', 'CHEBI:76579']
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
    'num_false_positives': 2,
    'num_true_negatives': 182401,
    'num_false_negatives': 31,
    'num_negatives': None,
    'precision': 0.9868421052631579,
    'recall': 0.8287292817679558,
    'f1': 0.9032258064516129,
    'accuracy': 0.9998182876992198
}
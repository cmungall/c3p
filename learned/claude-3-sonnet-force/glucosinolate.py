"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: CHEBI:16830 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    A glucosinolate has a central C atom bonded to a glycone group via S, a sulfonated oxime group via N, 
    and a side chain. The sulfate group and side chain have an anti stereochemical configuration across C=N.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glucose backbone pattern
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]([C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No glucose backbone found"
    
    # Look for sulfate group attached to glucose via sulfur
    sulfate_pattern = Chem.MolFromSmarts("[S]C(=N[O-,S]([O-])(=O)=O)")
    if not mol.HasSubstructMatch(sulfate_pattern, useChirality=True):
        return False, "Sulfate group not found or incorrect stereochemistry"
    
    # Look for side chain attached to central carbon
    side_chain_pattern = Chem.MolFromSmarts("[C]=[N]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if not side_chain_matches:
        return False, "No side chain found"
    
    # Check anti stereochemistry across C=N
    central_carbon = mol.GetAtomWithIdx(side_chain_matches[0][0])
    n_atom = mol.GetAtomWithIdx(side_chain_matches[0][1])
    if central_carbon.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW or n_atom.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Incorrect stereochemistry across C=N bond"
    
    return True, "Contains glucose backbone, sulfate group, and side chain with anti stereochemistry across C=N"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16830',
        'name': 'glucosinolate',
        'definition': 'Water-soluble anionic substituted thioglucosides. Glucosinolates have a central C atom which is bonded via an S atom to a glycone group and via an N atom to a sulfonated oxime group, and which also carries a side-group. The side-chain and sulfate group have an anti stereochemical configuration across the C=N double bond.',
        'parents': ['CHEBI:25537', 'CHEBI:18336']
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
    'num_true_positives': 237,
    'num_false_positives': 4,
    'num_true_negatives': 182405,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9834710743801653,
    'recall': 1.0,
    'f1': 0.9916796916796917,
    'accuracy': 0.9999782407631296
}
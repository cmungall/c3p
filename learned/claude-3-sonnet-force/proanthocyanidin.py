"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: CHEBI:18035 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavan-3-ol units (aromatic ring with OH at position 3 and heterocycle)
    flavan_pattern = Chem.MolFromSmarts("[cR1]1[cR1][cR1][cR1][cR1][cR2]1[C@@H]2[C@H]([O])c3c([O])c([O])c([O])c([O])c3[C@H]2[O]")
    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    
    # Look for oligomeric structure (two or more flavan units linked)
    if len(flavan_matches) < 2:
        return False, "Not oligomeric, less than 2 flavan-3-ol units found"
    
    # Check for ether or C-C linkages between units
    link_pattern = Chem.MolFromSmarts("[O][cR1][cR1]")
    link_matches = mol.GetSubstructMatches(link_pattern)
    
    if not link_matches:
        return False, "No ether or C-C linkages found between flavan units"
    
    # Count rotatable bonds to verify oligomeric nature
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for oligomer"
    
    # Check molecular weight - proanthocyanidins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for proanthocyanidin"
    
    return True, "Contains two or more hydroxyflavan units linked via ether or C-C bonds"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:18035',
        'name': 'proanthocyanidin',
        'definition': 'A flavonoid oligomer obtained by the the condensation of two or more units of hydroxyflavans.',
        'parents': ['CHEBI:16236', 'CHEBI:24851']
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
    'num_true_positives': 209,
    'num_false_positives': 4,
    'num_true_negatives': 182411,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.981308411214953,
    'recall': 1.0,
    'f1': 0.9905660377358491,
    'accuracy': 0.9998526141063597
}
"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:17794 flavone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavone(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is a flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavone skeleton pattern
    flavone_pattern = Chem.MolFromSmarts("c1cc(oc2ccccc2c1=O)-c")
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "Does not contain the flavone skeleton"
    
    # Count aromatic rings
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings != 3:
        return False, f"Expected 3 aromatic rings, found {aromatic_rings}"
    
    # Check for oxygens
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 4:
        return False, "Too few oxygens for a flavone"
    
    # Check for hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts("OC")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)
    if num_hydroxy < 1:
        return False, "No hydroxyl groups found"
    
    # Check for common flavone substituents
    substituents = ["CH3", "OCH3", "OH", "O", "CH2", "CH"]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            neighbors = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            substituent = "".join([get_symbol(n) for n in sorted(neighbors)])
            if substituent not in substituents:
                return False, f"Unusual substituent '{substituent}' found"
    
    return True, "Contains the flavone skeleton with expected aromatic rings, oxygens, and substituents"

def get_symbol(atomic_num):
    return Chem.Atom(atomic_num).GetSymbol()

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17794',
        'name': 'flavone',
        'definition': 'A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.',
        'parents': ['CHEBI:18594', 'CHEBI:35908']
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
    'num_true_positives': 151,
    'num_false_positives': 5,
    'num_true_negatives': 182412,
    'num_false_negatives': 16,
    'num_negatives': None,
    'precision': 0.9678379149631216,
    'recall': 0.9040404040404041,
    'f1': 0.9351303979560504,
    'accuracy': 0.9999530395199616
}
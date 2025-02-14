"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: CHEBI:49044 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from the triterpene hydrocarbon cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for tetracyclic scaffold
    tetracyclic_pattern = Chem.MolFromSmarts("[C&R1]1[C&R1]2[C&R1]3[C&R1]4[C&R1][C&R1][C&R1][C&R1]2[C&R1][C&R1][C&R1]4[C&R1][C&R1]3[C&R1]1")
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "Not a tetracyclic scaffold"
    
    # Check for cucurbitane core
    cucurbitane_pattern = Chem.MolFromSmarts("[C&R1]1[C&R1]2[C&R1]3[C&R1]4[C&R1][C&R1][C&R1][C&R1]2[C&R1][C&R1][C&R1]4[C&R1][C&R1]3[C&R1]1")
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "Does not contain a cucurbitane core"
    
    # Check for oxygenated substituents (hydroxy, carbonyl, ester, etc.)
    oxygenated_pattern = Chem.MolFromSmarts("[O;X2,X3]")
    oxygenated_matches = mol.GetSubstructMatches(oxygenated_pattern)
    if not oxygenated_matches:
        return False, "No oxygenated substituents found"
    
    # Check molecular weight - cucurbitacins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for cucurbitacin"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for cucurbitacin"
    if o_count < 3:
        return False, "Too few oxygens for cucurbitacin"
    
    return True, "Contains a tetracyclic cucurbitane core with oxygenated substituents"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:49044',
        'name': 'cucurbitacin',
        'definition': 'Any one of a class of tetracyclic triterpenoids, formally derived from the triterpene hydrocarbon cucurbitane, developed by some plants (especially those of the family Cucurbitaceaeas) as a defence mechanism against herbivores.',
        'parents': ['CHEBI:60027', 'CHEBI:63637']
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
    'num_true_positives': 56,
    'num_false_positives': 1,
    'num_true_negatives': 182414,
    'num_false_negatives': 13,
    'num_negatives': None,
    'precision': 0.9827586206896552,
    'recall': 0.8116883116883117,
    'f1': 0.8910891089108911,
    'accuracy': 0.9997794258063422
}
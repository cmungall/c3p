"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: CHEBI:27683 chlorophyll

A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings.
The rings can have various side chains which usually include a long phytol chain.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom
    has_mg = any(atom.GetAtomicNum() == 12 for atom in mol.GetAtoms())
    if not has_mg:
        return False, "No magnesium atom found"

    # Look for porphyrin ring pattern (4 pyrrole-like rings + 1 extra ring)
    porphyrin_pattern = Chem.MolFromSmarts("[N,n]!@[c,C]1[nH]c2[nH]c3[nH]c4[nH]c(c2c1c1c3c5c4c1c1c6c5c1c1c6c6cc6)c6cc6"
                                           "|!@[n,N]!@c1c2c3c4c5c6c1c2c3c4c5c6")
    
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin ring system found"

    # Look for long aliphatic side chain (phytol chain)
    phytol_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    phytol_matches = mol.GetSubstructMatches(phytol_pattern)
    if not phytol_matches:
        return False, "No long aliphatic side chain (phytol) found"

    # Check for other typical chlorophyll substituents (vinyl groups, esters, etc.)
    substituent_patterns = [
        Chem.MolFromSmarts("[C,c]=C=C"), # vinyl group
        Chem.MolFromSmarts("[CX3](=[OX1])OC"), # ester group
        Chem.MolFromSmarts("C(=O)C") # ketone
    ]
    has_substituents = any(mol.HasSubstructMatch(p) for p in substituent_patterns)
    if not has_substituents:
        return True, "Minimal chlorophyll structure: porphyrin ring system with Mg and aliphatic side chain"

    return True, "Contains porphyrin ring system with Mg, long aliphatic side chain, and typical chlorophyll substituents"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27683',
        'name': 'chlorophyll',
        'definition': 'A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.',
        'parents': ['CHEBI:26860', 'CHEBI:37083']
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
    'num_true_positives': 22,
    'num_false_positives': 0,
    'num_true_negatives': 1100,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}
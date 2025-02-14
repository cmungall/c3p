"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol bearing a single acyl, alkyl, or alk-1-enyl substituent
    at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern (three carbons with three hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Count the number of substituents on glycerol oxygen atoms
    # Get indices of oxygen atoms connected to carbons in glycerol
    glycerol_oxygens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 
                        and any(neigh.GetAtomicNum() == 6 for neigh in atom.GetNeighbors())]
    
    substituent_count = 0
    for o_idx in glycerol_oxygens:
        atom = mol.GetAtomWithIdx(o_idx)
        # Count non-hydrogen neighbors excluding the glycerol carbon
        neighbors = [neigh for neigh in atom.GetNeighbors() if neigh.GetAtomicNum() != 1]
        if len(neighbors) > 1:
            substituent_atom = [neigh for neigh in neighbors if neigh.GetAtomicNum() != 6][0]
            substituent_count += 1

    if substituent_count != 1:
        return False, f"Expected 1 substituent on glycerol, found {substituent_count}"

    # Check for acyl (ester linkage), alkyl (ether linkage), or alk-1-enyl (vinyl ether linkage)
    # Acyl group (ester linkage): glycerol oxygen connected to carbonyl carbon
    ester_pattern = Chem.MolFromSmarts("O[C;R0]=[O]")
    has_ester = mol.HasSubstructMatch(ester_pattern)

    # Alkyl group (ether linkage): glycerol oxygen connected to alkyl chain
    ether_pattern = Chem.MolFromSmarts("OCC")
    has_ether = mol.HasSubstructMatch(ether_pattern)

    # Alk-1-enyl group (vinyl ether linkage): glycerol oxygen connected to vinyl group
    vinyl_ether_pattern = Chem.MolFromSmarts("OC=C")
    has_vinyl_ether = mol.HasSubstructMatch(vinyl_ether_pattern)

    if has_ester:
        return True, "Contains glycerol backbone with a single acyl (ester-linked) substituent"
    elif has_ether:
        return True, "Contains glycerol backbone with a single alkyl (ether-linked) substituent"
    elif has_vinyl_ether:
        return True, "Contains glycerol backbone with a single alk-1-enyl (vinyl ether-linked) substituent"
    else:
        return False, "No valid substituent (acyl, alkyl, or alk-1-enyl) found on glycerol"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'monoradylglycerol',
        'definition': 'Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.',
        'parents': []
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
    'stdout': None
}
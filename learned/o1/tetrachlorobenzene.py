"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is any chlorobenzene carrying four chloro groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the ring information
    ri = mol.GetRingInfo()
    # Flag to indicate if tetrachlorobenzene ring is found
    found = False
    
    # Iterate over all atom rings
    for ring in ri.AtomRings():
        # Check if the ring is of size 6 (benzene ring)
        if len(ring) != 6:
            continue
        # Check if all atoms in the ring are aromatic carbons
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        is_carbon_ring = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)
        if not (is_aromatic and is_carbon_ring):
            continue
        
        # Count the number of chlorine substituents attached to the ring
        chlorine_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 17:
                    chlorine_count += 1
        if chlorine_count == 4:
            found = True
            break  # No need to check other rings
    
    if found:
        return True, "Contains benzene ring with exactly four chlorine substituents"
    else:
        return False, "Does not contain benzene ring with exactly four chlorine substituents"
    

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37384',
                              'name': 'tetrachlorobenzene',
                              'definition': 'Any member of the class of chlorobenzenes carrying four chloro groups at unspecified positions.',
                              'parents': ['CHEBI:33853']},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None,
        'num_true_positives': 150,
        'num_false_positives': 5,
        'num_true_negatives': 182400,
        'num_false_negatives': 20,
        'num_negatives': None,
        'precision': 0.967741935483871,
        'recall': 0.8823529411764706,
        'f1': 0.923076923076923,
        'accuracy': 0.999862189145545}
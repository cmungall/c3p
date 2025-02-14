"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: CHEBI:79067 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    A 3-substituted propionyl-CoA(4-) is an acyl-CoA(4-) where the acyl group is a propionyl
    group substituted at position 3 (beta carbon), attached via a thioester linkage to CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA core SMARTS pattern
    coa_smarts = "[#8]-[#6]-1-[#6]=[O]-[#6]=[N]-[#6]=[#7]-1"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA core structure not found"

    # Define thioester linkage pattern: sulfur connected to carbonyl carbon
    thioester_pattern = Chem.MolFromSmarts("S[C](=O)[C]")

    # Find matches for thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Assume first thioester match is the acyl linkage
    thioester_atoms = thioester_matches[0]
    sulfur_idx = thioester_atoms[0]
    carbonyl_idx = thioester_atoms[1]
    alpha_carbon_idx = thioester_atoms[2]

    # Get the beta carbon (3rd carbon from carbonyl carbon)
    # Traverse from alpha carbon excluding the carbonyl carbon
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    alpha_neighbors = [nbr.GetIdx() for nbr in alpha_carbon.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
    if not alpha_neighbors:
        return False, "No beta carbon found in acyl chain"

    beta_carbon_idx = alpha_neighbors[0]
    beta_carbon = mol.GetAtomWithIdx(beta_carbon_idx)

    # Check if beta carbon has substituents other than chain continuation and hydrogens
    beta_neighbors = beta_carbon.GetNeighbors()
    substituents = []
    for nbr in beta_neighbors:
        nbr_idx = nbr.GetIdx()
        if nbr_idx != alpha_carbon_idx and nbr.GetAtomicNum() != 1:
            substituents.append(nbr)

    if len(substituents) > 1:
        return True, "Beta carbon (3rd from carbonyl) is substituted"
    else:
        return False, "Beta carbon is unsubstituted"
    

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:79067',
                              'name': '3-substituted propionyl-CoA(4-)',
                              'definition': 'An acyl-CoA(4-) oxoanion arising from the '
                                            'deprotonation of the phosphate and diphosphate '
                                            'OH groups of any 3-substituted propionyl-CoA; '
                                            'major species at pH 7.3.',
                              'parents': ['CHEBI:57287', 'CHEBI:79066']},
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
        'num_false_positives': 4,
        'num_true_negatives': 182407,
        'num_false_negatives': 23,
        'num_negatives': None,
        'precision': 0.974025974025974,
        'recall': 0.8670520231213873,
        'f1': 0.9174311926605504,
        'accuracy': 0.9998521228585199}
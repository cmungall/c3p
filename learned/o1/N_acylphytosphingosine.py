"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:76955 N-acylphytosphingosine
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine
    having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amide bond pattern (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"
    
    # Define phytosphingosine backbone pattern
    # It consists of a chain with hydroxyl groups at positions 1,3,4 and connected to the amide nitrogen
    phytosphingosine_pattern = Chem.MolFromSmarts("""
    [N;X3]([C;X4;H2][O])[C;X4;H]([O])[C;X4;H]([O])[C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]
    """)
    
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Loop through all amide bonds to check for correct connections
    for match in amide_matches:
        carbonyl_c_idx = match[0]
        nitrogen_idx = match[2]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Check if nitrogen is connected to phytosphingosine backbone
        ns = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(nitrogen_idx, carbonyl_c_idx).GetIdx()], addDummies=False)))
        if ns.HasSubstructMatch(phytosphingosine_pattern):
            # Check acyl chain length (number of carbons in acyl group)
            acyl_chain = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(nitrogen_idx, carbonyl_c_idx).GetIdx()], addDummies=False)
            carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
            acyl_c_chain = [atom.GetIdx() for atom in Chem.rdmolops.GetShortestPath(mol, carbonyl_c_idx, carbonyl_c.GetNeighbors()[0].GetIdx()) if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6]
            acyl_chain_length = len(acyl_c_chain)
            if acyl_chain_length < 8:
                return False, "Acyl chain too short to be fatty acyl group"
            else:
                return True, "Molecule is an N-acylphytosphingosine"
    
    return False, "No matching phytosphingosine backbone attached via amide bond"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76955',
                              'name': 'N-acylphytosphingosine',
                              'definition': 'A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.',
                              'parents': ['CHEBI:75769']},
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
    'num_true_positives': 145,
    'num_false_positives': 5,
    'num_true_negatives': 182400,
    'num_false_negatives': 28,
    'num_negatives': None,
    'precision': 0.9668897816091954,
    'recall': 0.8389261744966443,
    'f1': 0.8987341772151899,
    'accuracy': 0.9998173668842193}
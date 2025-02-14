"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60940 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA is a fatty acyl-CoA that results from the formal condensation 
    of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for Coenzyme A (partial pattern to match key features)
    coA_smarts = Chem.MolFromSmarts("""
    NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO[P](=O)(O)O)[C@@H](O)[C@H]1O
    """)
    if not mol.HasSubstructMatch(coA_smarts):
        return False, "Coenzyme A moiety not found"

    # Define SMARTS pattern for thioester linkage between fatty acyl group and CoA
    thioester_smarts = Chem.MolFromSmarts("C(=O)SCCN")  # Simplified pattern including part of CoA
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Identify the acyl chain attached via thioester linkage
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon
        sulfur_idx = match[2]      # Sulfur atom

        # Get bond between carbonyl carbon and sulfur
        bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
        if bond is None:
            continue

        # Break the bond to isolate the fatty acyl chain
        mol_frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
        frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)
        
        # Identify the fragment containing the fatty acyl chain
        acyl_chain = None
        for frag in frags:
            atom_indices = [atom.GetIdx() for atom in frag.GetAtoms()]
            if carbonyl_c_idx in atom_indices:
                acyl_chain = frag
                break
        
        if acyl_chain is None:
            continue

        # Count the number of carbons in the acyl chain
        num_carbons = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)
        # Exclude the carbonyl carbon (assumed to be counted)
        num_carbons -= 1

        if num_carbons <= 5:
            return True, f"Found short-chain fatty acyl-CoA with acyl chain length {num_carbons}"
        else:
            return False, f"Acyl chain length is {num_carbons}, not short-chain"

    return False, "Thioester linkage not properly formed or acyl chain not identified"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:60940',
        'name': 'short-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.',
        'parents': ['CHEBI:37554', 'CHEBI:57288']
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
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}
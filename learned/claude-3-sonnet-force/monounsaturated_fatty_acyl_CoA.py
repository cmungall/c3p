"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:36691 monounsaturated fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is an unsaturated fatty acyl-CoA with one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA moiety
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Look for fatty acyl chain
    alkyl_pattern = Chem.MolFromSmarts("[CX4]([CX4H3])([CX4H3])([CX4H3])[CX4H2]")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No fatty acyl chain found"
    
    # Count double bonds
    num_double_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol) - sum(1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE)
    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds, expected 1"
    
    # Count carbon atoms in fatty acyl chain
    chain_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetIsAromatic() == 0]
    chain_length = sum(1 for atom in chain_atoms if atom.GetDegree() <= 2)
    if chain_length < 4:
        return False, "Fatty acyl chain too short"
    
    return True, "Contains a CoA moiety and a monounsaturated fatty acyl chain"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36691',
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.',
        'parents': ['CHEBI:35473', 'CHEBI:57737']
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
    'num_true_positives': 178,
    'num_false_positives': 8,
    'num_true_negatives': 182428,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.9573418025651016,
    'recall': 0.9675324675324675,
    'f1': 0.9624212654794595,
    'accuracy': 0.9998575796739161
}
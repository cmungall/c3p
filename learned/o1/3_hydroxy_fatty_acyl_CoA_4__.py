"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

    A 3-hydroxy fatty acyl-CoA(4-) is an acyl-CoA molecule where the fatty acyl chain has a hydroxyl group at the 3-position,
    and the molecule is deprotonated at the phosphate and diphosphate groups, resulting in a net charge of -4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check net charge
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is {total_charge}, expected -4"

    # Check for Coenzyme A moiety
    # Simplified CoA pattern
    coenzyme_a_smarts = 'NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@H]1OP([O-])([O-])=O'
    coenzyme_a = Chem.MolFromSmarts(coenzyme_a_smarts)
    if not mol.HasSubstructMatch(coenzyme_a):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage
    thioester_smarts = 'C(=O)SCCNC(=O)'
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "Thioester linkage not found"

    # Check for 3-hydroxy fatty acyl chain
    # The chain attached to the thioester should have a hydroxyl at the 3-position
    hydroxy_chain_smarts = '[C;!R][C@H](O)[C;!R][CX3](=O)SCCNC(=O)'
    hydroxy_chain = Chem.MolFromSmarts(hydroxy_chain_smarts)
    if not mol.HasSubstructMatch(hydroxy_chain):
        return False, "3-hydroxy fatty acyl chain not found"

    return True, "Molecule is a 3-hydroxy fatty acyl-CoA(4-)"

__metadata__ = {
    'chemical_class': {
        'name': '3-hydroxy fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.',
    },
    'success': True,
    'error': '',
}
"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

    A 3-hydroxy fatty acyl-CoA(4-) is an acyl-CoA molecule where the fatty acyl chain has a hydroxyl group at the 3-position,
    and the molecule is deprotonated at the phosphate and diphosphate OH groups, resulting in a net charge of -4.

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

    # Check for Coenzyme A moiety (more general SMARTS pattern)
    coenzyme_a_smarts = Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)O[P,O-]')
    if not mol.HasSubstructMatch(coenzyme_a_smarts):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage to fatty acyl chain
    thioester_smarts = Chem.MolFromSmarts('C(=O)SCCN')
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "Thioester linkage not found"

    # Check for 3-hydroxy fatty acyl chain attached to thioester
    # The chain attached to the thioester should have a hydroxyl at the 3-position from the carbonyl carbon
    hydroxy_acyl_chain_smarts = Chem.MolFromSmarts('C(=O)SC[C;!R][C@H](O)[C;!R]')
    if not mol.HasSubstructMatch(hydroxy_acyl_chain_smarts):
        return False, "3-hydroxy fatty acyl chain not found"

    # Check for deprotonated phosphate and diphosphate groups with negative charges
    phosphate_smarts = Chem.MolFromSmarts('P(=O)([O-])[O-]')
    num_phosphates = len(mol.GetSubstructMatches(phosphate_smarts))
    if num_phosphates < 3:
        return False, f"Found {num_phosphates} deprotonated phosphate groups, expected at least 3"

    # Ensure net charge is -4
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if total_charge != -4:
        return False, f"Total charge is {total_charge}, expected -4"

    return True, "Molecule is a 3-hydroxy fatty acyl-CoA(4-)"

__metadata__ = {
    'chemical_class': {
        'name': '3-hydroxy fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.',
    },
    'success': True,
    'error': '',
}
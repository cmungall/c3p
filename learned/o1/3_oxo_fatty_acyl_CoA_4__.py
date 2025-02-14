"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    A 3-oxo-fatty acyl-CoA(4-) is an acyl-CoA(4-) with a 3-oxo-fatty acyl chain attached via a thioester linkage,
    where the phosphate and diphosphate groups are deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for Coenzyme A (CoA) moiety
    # CoA SMARTS pattern (simplified version)
    coa_smarts = Chem.MolFromSmarts('NC(=O)CCNC(=O)C(O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](CN2C=NC3=C(N)N=CN=C32)[C@H](O)[C@@H]1OP(O)(O)=O')
    if not mol.HasSubstructMatch(coa_smarts):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage: C(=O)-S-
    thioester_smarts = Chem.MolFromSmarts('C(=O)SC')
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "Thioester linkage not found"

    # Check for beta-ketoacyl group
    # Pattern: O=C-C-C(=O)-C (beta-keto group in fatty acyl chain)
    beta_ketoacyl_smarts = Chem.MolFromSmarts('O=CCC(=O)C')
    if not mol.HasSubstructMatch(beta_ketoacyl_smarts):
        return False, "3-oxo (beta-ketoacyl) group not found in acyl chain"

    # Check for deprotonated phosphate groups
    # Count total negative charges (should be -4)
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Incorrect net charge ({total_charge}), expected -4 for deprotonated phosphate groups"

    # The molecule passes all checks
    return True, "Molecule is a 3-oxo-fatty acyl-CoA(4-) with deprotonated phosphate groups"

__metadata__ = {
    'chemical_class': {
        'name': '3-oxo-fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups of any 3-oxo-fatty acyl-CoA.',
        'examples': [
            '3-oxoisooctadecanoyl-CoA(4-)',
            '3-oxohexacosanoyl-CoA(4-)',
            '3-oxooctanoyl-CoA(4-)',
            # ... (additional examples)
        ]
    }
}
"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA_4__(smiles):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    A 3-oxo-fatty acyl-CoA(4-) is an acyl-CoA(4-) arising from deprotonation of the
    phosphate and diphosphate groups of any 3-oxo-fatty acyl-CoA.

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

    # Look for 3-oxo-fatty acid pattern
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)[C]")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo-fatty acid moiety found"

    # Look for CoA(4-) pattern
    coa_pattern = Chem.MolFromSmarts("C(C)(C)OP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA(4-) moiety found"

    # Look for thioester linkage between fatty acid and CoA
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found between fatty acid and CoA"

    # Count number of phosphate and diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[PX4]([O-])(=[O])([O-])")
    diphosphate_pattern = Chem.MolFromSmarts("[PX4]([O-])(=[O])([O-])O[PX4]([O-])(=[O])([O-])")
    n_phosphate = len(mol.GetSubstructMatches(phosphate_pattern))
    n_diphosphate = len(mol.GetSubstructMatches(diphosphate_pattern))
    if n_phosphate != 1 or n_diphosphate != 1:
        return False, "Incorrect number of phosphate/diphosphate groups"

    return True, "Contains 3-oxo-fatty acid moiety and CoA(4-) group linked via a thioester bond"
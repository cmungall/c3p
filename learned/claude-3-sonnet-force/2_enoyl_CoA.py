"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: CHEBI:36347 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA backbone"

    # Check for double bond between positions 2 and 3 of the acyl group
    enoyl_pattern = Chem.MolFromSmarts("C=CC(=O)SC")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No double bond between positions 2 and 3 of the acyl group"

    # Check for unsaturation in the acyl group
    acyl_group = Chem.MolFromSmiles(Chem.MolToSmiles(mol).split("SCCNC(=O)CCNC(=O)")[0] + "C(=O)O")
    if Chem.CalcNumRotatableBonds(acyl_group) == 0:
        return False, "The acyl group is saturated"

    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = acyl_group.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain in the acyl group"

    return True, "Contains an unsaturated fatty acyl-CoA with a double bond between positions 2 and 3"
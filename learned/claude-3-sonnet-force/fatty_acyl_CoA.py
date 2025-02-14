"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:33359 fatty acyl-CoA

A fatty acyl-CoA results from the formal condensation of the thiol group of coenzyme A 
with the carboxy group of any fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA substructure found"

    # Look for fatty acid substructure
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX3]")
    fatty_acid_match = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_match) == 0:
        return False, "No fatty acid substructure found"

    # Check for ester linkage between fatty acid and CoA
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[CX3]")
    ester_match = mol.GetSubstructMatches(ester_pattern)
    if len(ester_match) == 0:
        return False, "No ester linkage found between fatty acid and CoA"

    # Check for long aliphatic chain in fatty acid
    mcs = rdFMCS.FindMCS([mol, Chem.MolFromSmiles("CCCCCCCCCCCCCCC")], matchValences=True, completeRingsOnly=True)
    if mcs.numAtoms < 8:
        return False, "Fatty acid chain too short"

    return True, "Contains a fatty acid esterified to CoA"
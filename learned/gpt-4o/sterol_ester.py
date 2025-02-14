"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester linkage formed by condensation of a carboxylic acid
    with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for steroid skeleton and ester bond
    steroid_pattern = Chem.MolFromSmarts("[C;R1]1[C;R1][C;R1][C;R1]2[C;R1]3[C;R1][C;R1][C;R1]4[C;R1][C;R1][C;R1][C;R1]1[C;R1]2[C;R1]3[C;R1]4")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;R]")

    # Check for steroid skeleton
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid skeleton found"

    # Check for ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    return True, "Contains a steroid skeleton with an ester linkage"

# Test the function with a sterol ester example
smiles_example = "CCCCCCCCCCCCCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
is_sterol_ester(smiles_example)
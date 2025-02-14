"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: CHEBI:53685 tocol

A tocol is a chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a
saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chroman-6-ol skeleton
    chroman_pattern = Chem.MolFromSmarts("[OX2;R]1[CX4;R]2=[CX3]([CH3])[CX3]=[CX3][CX3]=[CX3]2[CX4]=[CX3][CX3]1")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No chroman-6-ol skeleton found"

    # Look for isoprenoid chain at position 2
    isoprenoid_pattern = Chem.MolFromSmarts("[CX4](C)(C)[CX4;R1]=[CX3][CX3]=[CX3][CX3]=[CX3][CX3]=[CX3][CX3]=[CX3]")
    if not mol.HasSubstructMatch(isoprenoid_pattern):
        return False, "No isoprenoid chain at position 2"

    # Count number of isoprenoid units
    isoprenoid_units = len(mol.GetSubstructMatches(isoprenoid_pattern))
    if isoprenoid_units != 3:
        return False, f"Found {isoprenoid_units} isoprenoid units, should be 3"

    # Check if chain is saturated or triply unsaturated
    unsaturated_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.TRIPLE)
    if unsaturated_bonds not in [0, 3]:
        return False, "Side chain should be either saturated or triply unsaturated"

    return True, "Contains a chroman-6-ol skeleton substituted at position 2 with a saturated or triply-unsaturated 3-isoprenoid chain"
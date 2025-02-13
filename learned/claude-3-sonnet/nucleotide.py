"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:25517 nucleotide

A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nucleobase pattern (aromatic rings with N atoms)
    base_pattern = Chem.MolFromSmarts("a1aaaanaa1")
    if not mol.HasSubstructMatch(base_pattern):
        return False, "No nucleobase found"

    # Look for sugar pattern (5-membered ring with 4 oxygen atoms)
    sugar_pattern = Chem.MolFromSmarts("OC1OC(O)C(O)C1O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Look for phosphate group (-O-P(=O)(O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycosidic bond between base and sugar
    glycosidic_bond = Chem.MolFromSmarts("[N,O]C1OC(O)C(O)C1O").GetSubstructMatch(mol)
    if not glycosidic_bond:
        return False, "No glycosidic bond between base and sugar"

    # Count phosphate groups (should be 1 or more)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"

    # Additional checks based on common nucleotide properties
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Too few rotatable bonds for nucleotide"

    n_rings = AllChem.CalcNumRings(mol)
    if n_rings < 2:
        return False, "Too few rings for nucleotide"

    n_hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    if n_hetero < 5:
        return False, "Too few heteroatoms for nucleotide"

    return True, "Contains nucleobase, sugar, and phosphate group"
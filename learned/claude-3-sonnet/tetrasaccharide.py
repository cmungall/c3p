"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:60598 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from typing import Tuple

def is_tetrasaccharide(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count monosaccharide units (assumes each monosaccharide has a ring)
    rings = mol.GetRingInfo().AtomRings()
    monosaccharide_units = sum(1 for ring in rings if len(ring) >= 5)
    if monosaccharide_units != 4:
        return False, f"Found {monosaccharide_units} monosaccharide units, need exactly 4"

    # Check for glycosidic bonds (O-C-O pattern)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_bond_matches) < 3:
        return False, f"Found {len(glycosidic_bond_matches)} glycosidic bonds, need at least 3"

    # Check molecular weight - tetrasaccharides typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for tetrasaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 16:
        return False, "Too few carbons for tetrasaccharide"
    if o_count < 12:
        return False, "Too few oxygens for tetrasaccharide"

    return True, "Contains four monosaccharide units linked by glycosidic bonds"
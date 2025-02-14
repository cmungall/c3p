"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
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

    # Count number of rings in molecule
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Check each ring for carbons and oxygens, count monosaccharide units
    monosaccharide_count = 0
    for ring in ring_info.AtomRings():
        carbon_count = 0
        oxygen_count = 0
        for atom_index in ring:
           atom = mol.GetAtomWithIdx(atom_index)
           if atom.GetAtomicNum() == 6:
               carbon_count += 1
           elif atom.GetAtomicNum() == 8:
               oxygen_count += 1
        if carbon_count >= 4 and oxygen_count > 0:
            monosaccharide_count += 1
    
    # Check for exactly 4 monosaccharide units
    if monosaccharide_count != 4:
        return False, f"Found {monosaccharide_count} monosaccharide units, requires 4 for tetrasaccharide"
    
    # Count glycosidic bonds (C-O-C connecting rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    # A tetrasaccharide will have at least 3 glycosidic bonds
    if len(glycosidic_bond_matches) < 3 :
         return False, f"Found {len(glycosidic_bond_matches)} glycosidic bonds, requires at least 3 for tetrasaccharide"

    return True, "Contains 4 monosaccharide units connected via glycosidic bonds"
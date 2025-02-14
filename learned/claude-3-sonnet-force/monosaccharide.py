"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:15955 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or polyhydroxy ketone with three or more carbon atoms,
    without glycosidic connection to other units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons, hydrogens, oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Must have at least 3 carbons
    if c_count < 3:
        return False, "Too few carbon atoms for monosaccharide"

    # Check for carbonyl group (aldehyde or ketone)
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found, not an aldehyde or ketone"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Check for polyhydroxy (at least 3 hydroxyls)
    if hydroxyl_count < 3:
        return False, "Not enough hydroxyl groups for monosaccharide"

    # Check for glycosidic bonds (shouldn't have any)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Contains glycosidic bonds, not a monosaccharide"

    # Check hydrogen deficiency (for ring structures)
    hd = rdMolDescriptors.CalcHydrogenDeficiency(mol)
    if hd > 2:
        return False, "Hydrogen deficiency too high, likely polycyclic"

    # Additional checks based on counts
    if h_count < 2 * o_count:
        return False, "Too few hydrogens for the number of oxygens"
    if c_count > 8:
        return False, "Carbon chain too long for monosaccharide"

    return True, "Contains carbonyl group and at least 3 hydroxyls without glycosidic bonds"
"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: monosaccharide
Determines if a molecule is a monosaccharide based on its SMILES string.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with three or more carbon atoms,
    and without glycosidic connections to other saccharide units. It includes cyclic hemiacetals
    and hemiketals (furanose and pyranose forms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check number of carbon atoms (must be 3 or more)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Contains {c_count} carbon atoms, which is less than 3"

    # Check for aldehyde group (terminal carbonyl group)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Check for ketone group (internal carbonyl group)
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # Check for cyclic hemiacetal or hemiketal (ring forms)
    # Five or six-membered ring with one oxygen (furanose or pyranose)
    ring_info = mol.GetRingInfo()
    ring_matches = []
    for ring in ring_info.AtomRings():
        if len(ring) == 5 or len(ring) == 6:
            oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygens_in_ring == 1:
                ring_matches.append(ring)
    has_cyclic_hemiacetal = len(ring_matches) > 0

    if not (has_aldehyde or has_ketone or has_cyclic_hemiacetal):
        return False, "Does not contain an aldehyde, ketone, or cyclic hemiacetal/hemiketal group"

    # Check for multiple hydroxyl groups (-OH) attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Contains {len(hydroxyl_matches)} hydroxyl groups attached to carbons, which is less than 2"

    # Ensure molecule is a single saccharide unit (no glycosidic bonds)
    # Glycosidic bonds are acetal linkages at anomeric carbon connecting two saccharides
    # Anomeric carbon is acetal (connected to two oxygens) in glycosidic bonds
    acetal_pattern = Chem.MolFromSmarts("[C;D2](-[O;H1])[O;!H1]")
    if mol.HasSubstructMatch(acetal_pattern):
        return False, "Contains glycosidic bonds indicating connection to other saccharide units"

    # Count number of rings (should be 0 or 1 for monosaccharides)
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return False, f"Contains {num_rings} rings, which is more than expected for a monosaccharide"

    # Check that molecule does not have unusual elements (only C, H, O, N, S, P)
    allowed_atomic_nums = {1, 6, 7, 8, 15, 16}  # H, C, N, O, P, S
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains element with atomic number {atom.GetAtomicNum()}, which is not typical in monosaccharides"

    # Optional: Exclude very large molecules based on heavy atom count
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 25:
        return False, f"Contains {heavy_atom_count} heavy atoms, which is too large for a monosaccharide"

    return True, "Molecule satisfies criteria for a monosaccharide"
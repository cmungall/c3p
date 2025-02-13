"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon
    atoms, typically with multiple hydroxyl groups and potentially in cyclic form.

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

    # Check if the molecule has at least 3 carbon atoms, which is a minimum for monosaccharides
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms for a monosaccharide"

    # Consider the structure as potentially cyclic, check for a common ring size like 5 (furanose) or 6 (pyranose)
    rings = mol.GetRingInfo().AtomRings()
    has_valid_ring = any(len(ring) in {5, 6} for ring in rings)

    # Check for oxygens to infer hydroxyl groups presence in sufficient amount
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < min(c_count / 2, 3):  # at least half the number of carbons are hydroxylated
        return False, "Too few oxygen atoms for typical monosaccharide hydroxylation"

    # Look for a generic carbonyl group pattern or infer potential from ring structure
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    if not mol.HasSubstructMatch(carbonyl_pattern) and not has_valid_ring:
        return False, "Does not contain a carbonyl group or valid cyclic form"

    # Avoid overly complex molecules that are likely not a single monosaccharide molecule
    if c_count > 10 or o_count > 6:  # Arbitrary threshold for simplicity
        return False, "Too complex structures have been ruled out for clarity"

    return True, "Structure has characteristic features of monosaccharides (polyhydroxy and cyclic or carbonyl)"
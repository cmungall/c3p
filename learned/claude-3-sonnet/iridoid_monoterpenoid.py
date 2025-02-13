"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:25610 iridoid monoterpenoids

Iridoid monoterpenoids are monoterpenoids biosynthesized from isoprene and often intermediates
in the biosynthesis of alkaloids. They typically consist of a cyclopentane ring fused to a
six-membered oxygen heterocycle. Cleavage of a bond in the cyclopentane ring gives rise to
the subclass known as secoiridoids.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclopentane-oxygen heterocycle pattern
    iridoid_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2[C@@H]([C@H](C[C@@]1(C)O2)C)C")
    if not mol.HasSubstructMatch(iridoid_pattern):
        return False, "No iridoid core structure found"

    # Check if monoterpenoid (C10 skeleton)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons != 10:
        return False, "Not a monoterpenoid (skeleton does not contain 10 carbons)"

    # Check for rings and ring sizes
    rings = mol.GetRingInfo().AtomRings()
    ring_sizes = [len(ring) for ring in rings]
    if 5 not in ring_sizes or 6 not in ring_sizes:
        return False, "Lacks cyclopentane and/or 6-membered oxygen heterocycle"

    # Check for oxygen atoms
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 2:
        return False, "Fewer than 2 oxygens (iridoids typically have 2-4)"

    # Todo: Add more specific checks for common iridoid substituents/decorations

    return True, "Contains iridoid core structure (cyclopentane fused to 6-membered oxygen heterocycle)"
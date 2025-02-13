"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: CHEBI:33738 cyclohexenone
A cyclohexenone is any six-membered alicyclic ketone having one double bond in the ring.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenone(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one ring with 6 atoms
    ring_info = mol.GetRingInfo()
    rings = [ring for ring in ring_info.AtomRings() if len(ring) == 6]
    if len(rings) != 1:
        return False, "Must have exactly one 6-membered ring"

    # Check if the ring is alicyclic (no heteroatoms in the ring)
    ring_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[0]]
    if any(atom.GetAtomicNum() != 6 for atom in ring_atoms):
        return False, "Ring must be alicyclic (no heteroatoms)"

    # Check for exactly one C=O group in the ring
    carbonyl_pattern = Chem.MolFromSmarts("[C=O]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    carbonyl_ring_matches = [match for match in carbonyl_matches if mol.GetAtomWithIdx(match) in ring_atoms]
    if len(carbonyl_ring_matches) != 1:
        return False, "Must have exactly one carbonyl group in the ring"

    # Check for exactly one C=C bond in the ring
    double_bond_pattern = Chem.MolFromSmarts("[C=C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    double_bond_ring_matches = [match for match in double_bond_matches if mol.GetBondBetweenAtoms(match[0], match[1]).IsInRing()]
    if len(double_bond_ring_matches) != 1:
        return False, "Must have exactly one double bond in the ring"

    return True, "Molecule is a cyclohexenone"
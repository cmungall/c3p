"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as compounds having a fully conjugated cyclic dione structure,
    derived from aromatic compounds by conversion of an even number of -CH= groups into
    -C(=O)- groups with any necessary rearrangement of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        return False, "Molecule does not contain any rings"

    # Find all carbonyl carbons (C=O)
    carbonyl_carbons = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or \
               (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6):
                if begin_atom.GetAtomicNum() == 6:
                    carbon_atom = begin_atom
                else:
                    carbon_atom = end_atom
                carbonyl_carbons.add(carbon_atom.GetIdx())

    # Check each ring for quinone pattern
    for ring in atom_rings:
        ring_set = set(ring)
        # Check if ring is fully conjugated (all atoms are sp2 hybridized)
        is_conjugated = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.HybridizationType.SP2:
                is_conjugated = False
                break
        if not is_conjugated:
            continue  # Skip this ring

        # Count carbonyl carbons in ring
        carbonyls_in_ring = ring_set.intersection(carbonyl_carbons)
        num_carbonyls_in_ring = len(carbonyls_in_ring)
        if num_carbonyls_in_ring >= 2:
            return True, "Molecule contains a fully conjugated cyclic diketone characteristic of quinones"

    return False, "Molecule does not contain the characteristic quinone structure"
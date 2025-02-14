"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is any six-membered alicyclic ketone having one double bond in the ring.

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

    # Get ring info
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate cyclohexenone found
    cyclohexenone_found = False

    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if ring is six-membered
        if len(ring) != 6:
            continue

        # Check if ring is alicyclic (non-aromatic)
        is_alicyclic = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                is_alicyclic = False
                break
        if not is_alicyclic:
            continue

        # Collect bonds in the ring
        bond_idxs = []
        for i in range(len(ring)):
            atom1_idx = ring[i]
            atom2_idx = ring[(i + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if bond is not None:
                bond_idxs.append(bond.GetIdx())

        # Count double bonds in the ring
        double_bond_count = 0
        for bond_idx in bond_idxs:
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

        if double_bond_count != 1:
            continue  # Need exactly one double bond in the ring

        # Check for ketone group in the ring
        ketone_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Skip if not carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue  # Skip if neighbor is in the ring
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone_found = True
                        break
            if ketone_found:
                break

        if ketone_found:
            cyclohexenone_found = True
            break  # No need to check other rings

    if cyclohexenone_found:
        return True, "Contains a six-membered alicyclic ring with one double bond and a ketone group"
    else:
        return False, "Does not contain a cyclohexenone ring structure"
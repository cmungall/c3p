"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: Catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component,
    which is an aromatic ring with two hydroxyl groups in adjacent positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains catechol moiety, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()

    # Iterate over all rings in the molecule
    for ring in ri.AtomRings():
        # Check if the ring is aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Create a set for quick lookup of ring atom indices
            ring_atoms = set(ring)
            # Find atoms in the ring with attached hydroxyl groups
            hydroxylated_atoms = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Check if the atom is a carbon with degree >= 2
                if atom.GetAtomicNum() == 6 and atom.GetDegree() >= 2:
                    for neighbor in atom.GetNeighbors():
                        # Check for hydroxyl group (oxygen with single bond and degree 1)
                        if (neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1 and
                            mol.GetBondBetweenAtoms(idx, neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE):
                            hydroxylated_atoms.append(idx)
            # Check for adjacent hydroxylated carbons in the ring
            ring_size = len(ring)
            for i in range(ring_size):
                idx1 = ring[i]
                idx2 = ring[(i + 1) % ring_size]  # Next atom in the ring (with wrap-around)
                if idx1 in hydroxylated_atoms and idx2 in hydroxylated_atoms:
                    return True, "Contains catechol moiety (o-diphenol component)"

    return False, "Does not contain catechol moiety (o-diphenol component)"
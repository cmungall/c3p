"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:25122 pyrroline
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is any organic heteromonocyclic compound with a structure based on dihydropyrrole,
    which is a five-membered, non-aromatic ring containing one nitrogen atom and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    # Get list of atom indices for each ring
    rings = ring_info.AtomRings()

    # Initialize flag
    is_pyrroline = False

    # Iterate over each ring
    for ring in rings:
        # Check if ring is five-membered
        if len(ring) != 5:
            continue

        # Collect atoms in the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Count the number of nitrogen atoms in the ring
        num_nitrogens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        if num_nitrogens != 1:
            continue

        # Check if ring is aromatic
        is_aromatic = any(atom.GetIsAromatic() for atom in ring_atoms)
        if is_aromatic:
            continue

        # Count the number of double bonds in the ring
        num_double_bonds = 0
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in ring and end_idx in ring:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    num_double_bonds += 1

        if num_double_bonds != 1:
            continue

        # If all conditions are met, it's a pyrroline
        is_pyrroline = True
        break

    if is_pyrroline:
        return True, "Contains a non-aromatic five-membered ring with one nitrogen atom and one double bond (pyrroline ring)"
    else:
        return False, "Does not contain a pyrroline ring"
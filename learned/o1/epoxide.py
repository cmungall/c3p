"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:32987 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Iterate over each ring in the molecule
    for ring in atom_rings:
        # Check if the ring has exactly 3 atoms
        if len(ring) == 3:
            oxygen_count = 0
            carbon_count = 0
            # Check each atom in the ring
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                atomic_num = atom.GetAtomicNum()
                if atomic_num == 8:
                    # Atom is oxygen
                    oxygen_count += 1
                elif atomic_num == 6:
                    # Atom is carbon
                    carbon_count += 1
                else:
                    # Atom is neither carbon nor oxygen, not an epoxide ring
                    break
            else:
                # All atoms in ring have been checked
                if oxygen_count == 1 and carbon_count == 2:
                    return True, "Contains an epoxide ring (3-membered ring with one oxygen atom)"
    # No epoxide ring found
    return False, "Does not contain an epoxide ring"
"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl (PCB) based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings, with no other substituents or functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that molecule contains only C, H, and Cl atoms
    allowed_atomic_nums = {6, 1, 17}  # Carbon, Hydrogen, Chlorine
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom other than C, H, or Cl: {atom.GetSymbol()}"

    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_atom_sets = ring_info.AtomRings()

    # Find aromatic rings
    aromatic_rings = []
    for ring in ring_atom_sets:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))

    # Look for biphenyl core (two connected benzene rings)
    found_biphenyl = False
    for i in range(len(aromatic_rings)):
        for j in range(i + 1, len(aromatic_rings)):
            ring1 = aromatic_rings[i]
            ring2 = aromatic_rings[j]
            # Check if rings share a bond (are connected)
            connected = False
            for atom_idx1 in ring1:
                for atom_idx2 in ring2:
                    bond = mol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
                    if bond is not None:
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            connected = True
                            break
                if connected:
                    break
            if connected:
                # Found biphenyl core
                biphenyl_atoms = ring1.union(ring2)
                # Check substituents on biphenyl carbons
                chloro_count = 0
                is_valid = True
                for atom_idx in biphenyl_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx not in biphenyl_atoms:
                            if neighbor.GetAtomicNum() == 17:  # Chlorine
                                chloro_count += 1
                            else:
                                is_valid = False
                                break
                    if not is_valid:
                        break
                if not is_valid:
                    continue  # Invalid biphenyl, check next possibility
                # Check that all atoms are part of biphenyl core or are chlorines attached to it
                biphenyl_and_cl_atoms = biphenyl_atoms.copy()
                for atom_idx in biphenyl_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx not in biphenyl_atoms:
                            biphenyl_and_cl_atoms.add(neighbor_idx)
                if len(biphenyl_and_cl_atoms) != mol.GetNumAtoms():
                    continue  # Molecule contains other atoms not part of biphenyl or chlorines
                # Check chlorine count
                if 2 <= chloro_count <= 10:
                    return True, f"Contains biphenyl core with {chloro_count} chlorine atoms attached"
                else:
                    return False, f"Contains biphenyl core with {chloro_count} chlorine atoms, which is not between 2 and 10"
    return False, "No valid polychlorobiphenyl core found"
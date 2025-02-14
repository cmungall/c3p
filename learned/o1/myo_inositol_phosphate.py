"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify all cyclohexane rings (potential inositol rings)
    ssr = Chem.GetSymmSSSR(mol)
    cyclohexane_rings = []
    for ring in ssr:
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetAtomicNum() == 6 for atom in ring_atoms):
                cyclohexane_rings.append(ring)

    if not cyclohexane_rings:
        return False, "No cyclohexane rings found"

    found_inositol_phosphate = False
    for ring in cyclohexane_rings:
        is_inositol = True
        ring_atom_indices = set(ring)
        connected_atoms = set(ring)
        phosphate_found = False

        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)

            # Each carbon in the ring should have exactly three bonds:
            # two to other carbons in the ring and one to oxygen
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 3:
                is_inositol = False
                break

            # Identify the oxygen attached to this carbon
            oxygen_atom = None
            for neighbor in neighbors:
                n_idx = neighbor.GetIdx()
                if neighbor.GetAtomicNum() == 8 and n_idx not in ring_atom_indices:
                    oxygen_atom = neighbor
                    connected_atoms.add(n_idx)
                    break
            if oxygen_atom is None:
                is_inositol = False
                break

            # Check that oxygen is connected either to hydrogen or phosphorus
            oxygen_neighbors = oxygen_atom.GetNeighbors()
            non_carbon_neighbors = [n for n in oxygen_neighbors if n.GetIdx() != idx]
            if len(non_carbon_neighbors) > 1:
                is_inositol = False
                break
            elif len(non_carbon_neighbors) == 1:
                # Oxygen connected to phosphorus (phosphate group)
                phosphorus_atom = non_carbon_neighbors[0]
                if phosphorus_atom.GetAtomicNum() != 15:
                    is_inositol = False
                    break
                connected_atoms.add(phosphorus_atom.GetIdx())
                # Check that phosphorus is part of phosphate group
                p_neighbors = phosphorus_atom.GetNeighbors()
                phosphate_oxygens = [pn for pn in p_neighbors if pn.GetAtomicNum() == 8 and pn.GetIdx() != oxygen_atom.GetIdx()]
                if len(phosphate_oxygens) < 3:
                    is_inositol = False
                    break
                connected_atoms.update([pn.GetIdx() for pn in phosphate_oxygens])
                phosphate_found = True
            # If oxygen is only connected to carbon and hydrogen, that's acceptable (OH group)

        if is_inositol:
            # Ensure that all atoms in the molecule are connected to the ring via the oxygens
            # Collect all atoms connected to the inositol ring
            mol_atom_indices = set(range(mol.GetNumAtoms()))
            if connected_atoms != mol_atom_indices:
                is_inositol = False
                break
            # If we reach here, we have found a myo-inositol phosphate
            found_inositol_phosphate = True
            break

    if not found_inositol_phosphate:
        return False, "No myo-inositol phosphate structure found"

    return True, "Contains myo-inositol phosphate core structure"
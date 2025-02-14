"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES with special handling for pentavalent phosphorus
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string"

    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
    except:
        return False, "Invalid SMILES string"

    # Identify phosphorus atoms
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atoms found"

    terminal_phosphate_found = False
    for p_atom in p_atoms:
        # Check if P atom is part of a phosphate or diphosphate group
        o_neighbors = [nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(o_neighbors) < 3:
            continue  # Not a phosphate group
        # Check for double bonded oxygen
        num_pdbonds = sum(1 for bond in p_atom.GetBonds()
                          if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and
                          bond.GetOtherAtom(p_atom).GetAtomicNum() == 8)
        if num_pdbonds != 1:
            continue  # Not a phosphate group
        # Look for O atom connected to carbon
        for o_atom in o_neighbors:
            if o_atom.GetDegree() == 2:
                # O atom connected to P and possibly C
                for o_nbr in o_atom.GetNeighbors():
                    if o_nbr.GetIdx() == p_atom.GetIdx():
                        continue
                    if o_nbr.GetAtomicNum() == 6:
                        # Found phosphate attached via O to carbon
                        c_atom = o_nbr
                        terminal_phosphate_found = True
                        break
                if terminal_phosphate_found:
                    break
        if terminal_phosphate_found:
            break

    if not terminal_phosphate_found:
        return False, "Phosphate group not attached via oxygen to carbon"

    # Get the fragment (chain) attached to the carbon atom
    def get_fragment(mol, atom_idx):
        fragment = set()
        atoms_to_visit = [atom_idx]
        while atoms_to_visit:
            idx = atoms_to_visit.pop()
            if idx not in fragment:
                fragment.add(idx)
                atom = mol.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if neighbor.GetAtomicNum() == 6 and n_idx not in fragment:
                        atoms_to_visit.append(n_idx)
        return fragment

    fragment = get_fragment(mol, c_atom.GetIdx())

    num_carbons = len(fragment)
    if num_carbons < 15:
        return False, f"Chain connected to phosphate group is too short: {num_carbons} carbons"

    # Count double bonds in the fragment
    double_bonds = 0
    for idx in fragment:
        atom = mol.GetAtomWithIdx(idx)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if begin_idx in fragment and end_idx in fragment:
                    double_bonds += 1

    if double_bonds < 2:
        return False, f"Not enough double bonds in the chain: {double_bonds}"

    return True, "Molecule is a polyprenol phosphate with a polyprenol backbone and terminal phosphate group attached via oxygen"
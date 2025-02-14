"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a long-chain isoprenoid alcohol composed of more than one isoprene unit.
    They have the general formula H-[CH2C(Me)=CHCH2]nOH, with n > 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal primary alcohol group (O-H connected to a primary carbon)
    alcohol_pattern = Chem.MolFromSmarts('[CX4][OX2H]')
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No terminal primary alcohol group found"

    # Define a more general isoprene unit pattern
    isoprene_pattern = Chem.MolFromSmarts('C(=C)CC[C]')
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern, uniquify=True)
    num_isoprene_units = len(isoprene_matches)

    if num_isoprene_units < 2:
        return False, f"Found {num_isoprene_units} isoprene units, need more than one"

    # Check that isoprene units are connected in a chain (head-to-tail)
    # Create a set of atoms involved in isoprene units
    isoprene_atoms = set()
    for match in isoprene_matches:
        isoprene_atoms.update(match)

    # Check connectivity
    visited_atoms = set()
    queue = [alcohol_matches[0][0]]  # Start from the terminal alcohol carbon
    while queue:
        atom_idx = queue.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in isoprene_atoms and neighbor_idx not in visited_atoms:
                queue.append(neighbor_idx)

    # Ensure all isoprene atoms are connected
    if not isoprene_atoms.issubset(visited_atoms):
        return False, "Isoprene units are not connected in a chain"

    return True, f"Contains terminal alcohol and {num_isoprene_units} isoprene units connected in a chain"
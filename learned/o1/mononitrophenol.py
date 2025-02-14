"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    aromatic_rings = [ring for ring in atom_rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]

    # Check that there is exactly one aromatic ring
    if len(aromatic_rings) != 1:
        return False, f"Molecule has {len(aromatic_rings)} aromatic rings, expected 1"

    ring_atoms = set(aromatic_rings[0])

    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)

    # Check for exactly one nitro group
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, expected 1"

    # Check that nitro group is attached to the aromatic ring
    nitro_nitrogen_idx = nitro_matches[0][0]
    nitro_nitrogen_atom = mol.GetAtomWithIdx(nitro_nitrogen_idx)
    nitro_connected = False
    for neighbor in nitro_nitrogen_atom.GetNeighbors():
        if neighbor.GetIdx() in ring_atoms:
            nitro_connected = True
            break
    if not nitro_connected:
        return False, "Nitro group is not attached to the aromatic ring"

    # Define phenolic hydroxyl group pattern
    phenol_pattern = Chem.MolFromSmarts('[c][O;H1]')
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)

    # Check for exactly one phenolic hydroxyl group
    if len(phenol_matches) != 1:
        return False, f"Found {len(phenol_matches)} phenolic OH groups, expected 1"

    # Check that hydroxyl group is attached to the aromatic ring
    phenol_oxygen_idx = phenol_matches[0][1]
    phenol_oxygen_atom = mol.GetAtomWithIdx(phenol_oxygen_idx)
    hydroxyl_connected = False
    for neighbor in phenol_oxygen_atom.GetNeighbors():
        if neighbor.GetIdx() in ring_atoms:
            hydroxyl_connected = True
            break
    if not hydroxyl_connected:
        return False, "Phenolic OH group is not attached to the aromatic ring"

    # Check for other substituents on the aromatic ring
    allowed_substituents = {'C', 'H', 'O', 'N', 'F', 'Cl', 'Br', 'I'}
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in ring_atoms:
                element = neighbor.GetSymbol()
                if element not in allowed_substituents:
                    return False, f"Disallowed substituent '{element}' attached to aromatic ring"
                # Exclude large substituents (more than 3 atoms)
                substituent = Chem.PathToSubmol(mol, [idx, neighbor_idx])
                if substituent.GetNumHeavyAtoms() > 4:
                    return False, "Substituent attached to ring is too large"

    # Check that there are no other rings in the molecule
    if len(atom_rings) != 1:
        return False, f"Molecule has {len(atom_rings)} rings, expected 1"

    return True, "Molecule is a mononitrophenol with single aromatic ring and appropriate substituents"
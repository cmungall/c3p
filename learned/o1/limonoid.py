"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a triterpenoid that is highly oxygenated and has a prototypical structure
    containing or derived from a 4,4,8-trimethyl-17-furanylsteroid skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Check for at least four rings
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    if num_rings < 4:
        return False, f"Only {num_rings} rings found, need at least 4 rings"

    # Step 2: Check for furan ring
    furan_pattern = Chem.MolFromSmarts('c1ccoc1')
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found in molecule"

    # Step 3: Check if furan ring is connected to the ring system
    ring_atom_indices = set()
    for ring in ri.AtomRings():
        ring_atom_indices.update(ring)
    furan_atom_indices = set()
    for match in furan_matches:
        furan_atom_indices.update(match)
    # Find if any atom in the furan ring is connected to any atom in the other rings
    connected = False
    for idx in furan_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in ring_atom_indices and neighbor.GetIdx() not in furan_atom_indices:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Furan ring is not connected to the core ring system"

    # Step 4: Check for high oxygenation (e.g., at least 6 oxygen atoms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Not highly oxygenated, only {o_count} oxygen atoms found"

    # Step 5: Check for multiple methyl groups attached to ring carbons
    # Count methyl groups attached to ring carbons
    methyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.IsInRing():
                methyl_count += 1
    if methyl_count < 3:
        return False, f"Only {methyl_count} methyl groups attached to ring carbons, need at least 3"

    # Step 6: Check for triterpenoid skeleton by counting carbon atoms
    # Triterpenoids typically have around 30 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:
        return False, f"Only {c_count} carbon atoms found, likely not a triterpenoid"

    return True, "Molecule is classified as a limonoid"
"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.

    A 3'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring (B ring) attached at position 2 of the flavanone core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS pattern (2-phenylchroman-4-one)
    flavanone_core_smarts = 'O=C1CCC2=C(C1)C=CC=C2-c1ccccc1'
    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Flavanone core not found"

    # Get matches for the flavanone core
    core_matches = mol.GetSubstructMatches(flavanone_core)
    if not core_matches:
        return False, "Flavanone core not found"

    # Assume first match is the flavanone core
    core_atoms = core_matches[0]

    # Atom indices in the flavanone core SMARTS (0-based)
    # [O]=C1CCC2=C(C1)C=CC=C2-c1ccccc1
    # Index 7 is the attachment point to the B ring (phenyl ring)
    attachment_idx_in_smarts = 7
    attachment_atom_idx = core_atoms[attachment_idx_in_smarts]

    # Get the B ring (phenyl ring attached at position 2)
    attachment_atom = mol.GetAtomWithIdx(attachment_atom_idx)
    b_ring = None
    for neighbor in attachment_atom.GetNeighbors():
        if neighbor.GetIdx() not in core_atoms:
            # This neighbor is part of the B ring
            if neighbor.IsInRing():
                b_ring = neighbor.GetOwningMol().GetRingInfo().AtomRings()
                break
    if b_ring is None:
        return False, "B ring not found"

    # Find the B ring containing the neighbor atom
    for ring in b_ring:
        if neighbor.GetIdx() in ring:
            b_ring_atoms = ring
            break
    else:
        return False, "B ring not found"

    # Reorder B ring atoms starting from the attachment point
    b_ring_atoms = list(b_ring_atoms)
    attachment_idx_in_ring = b_ring_atoms.index(neighbor.GetIdx())
    b_ring_atoms = b_ring_atoms[attachment_idx_in_ring:] + b_ring_atoms[:attachment_idx_in_ring]

    # Position 1' is the attachment point, positions increase around the ring
    # Position 3' is at index 2
    pos_3_prime_idx = b_ring_atoms[2]
    pos_3_prime_atom = mol.GetAtomWithIdx(pos_3_prime_idx)

    # Check for hydroxy group at position 3'
    hydroxy_found = False
    for nbr in pos_3_prime_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
            hydroxy_found = True
            break

    if hydroxy_found:
        return True, "Contains flavanone core with hydroxy group at 3' position of B ring"
    else:
        return False, "No hydroxy group at 3' position of B ring"

__metadata__ = {
    'chemical_class': {
        'name': "3'-hydroxyflavanones",
        'definition': "Any hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring."
    }
}
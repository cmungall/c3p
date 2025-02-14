"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:35480 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde group (-C=O)
    aldehyde_pattern = Chem.MolFromSmarts("[C;H1]=[O;H0]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check if aldehyde is attached to an aromatic ring
    for match in aldehyde_matches:
        aromatic_atom_idx = None
        for atom_idx in mol.GetAtomWithIdx(match).GetNeighbors():
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsAromatic():
                aromatic_atom_idx = atom_idx
                break

        if aromatic_atom_idx is not None:
            # Get aromatic ring info
            ring_info = mol.GetRingInfo().AtomRings()
            for ring in ring_info:
                if aromatic_atom_idx in ring:
                    return True, "Aldehyde group attached to an aromatic ring"

    return False, "Aldehyde group not attached to an aromatic ring"
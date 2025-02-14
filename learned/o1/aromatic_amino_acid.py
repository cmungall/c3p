"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is an amino acid whose structure includes an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify amino acid backbone: [N][C][C(=O)O]
    amino_acid_pattern = Chem.MolFromSmarts("[N;D2]-[C;D3]-[C](=O)-O")
    aa_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not aa_matches:
        return False, "No amino acid backbone found"

    # Check for aromatic ring
    aromatic_rings = mol.GetAromaticRings()
    if not aromatic_rings:
        return False, "No aromatic rings found in molecule"

    # Check if aromatic ring is part of the side chain connected to alpha carbon
    for match in aa_matches:
        n_idx, alpha_c_idx, c_idx = match  # Indices of N, alpha C, and carbonyl C

        # Find side chain attached to alpha carbon
        alpha_c = mol.GetAtomWithIdx(alpha_c_idx)
        side_chain_atoms = []
        for neighbor in alpha_c.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx != n_idx and nbr_idx != c_idx:
                # Traverse side chain to find if it contains aromatic atoms
                visited = set()
                stack = [nbr_idx]
                while stack:
                    atom_idx = stack.pop()
                    if atom_idx in visited:
                        continue
                    visited.add(atom_idx)
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetIsAromatic():
                        return True, "Aromatic ring found in side chain"
                    for nbr in atom.GetNeighbors():
                        nbr_idx2 = nbr.GetIdx()
                        if nbr_idx2 != alpha_c_idx and nbr_idx2 not in visited:
                            stack.append(nbr_idx2)

    return False, "Aromatic ring not part of side chain connected to amino acid backbone"

__metadata__ = {   'chemical_class': {   'name': 'aromatic amino acid',
                                         'definition': 'An amino acid whose structure includes an aromatic ring.'}}
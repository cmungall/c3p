"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI:35457 monoamine

A monoamine is an aralylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.
Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones by the action of aromatic amino acid decarboxylase enzymes.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a monoamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aromatic ring
    aromatic_ring = mol.GetAromaticRings()
    if not aromatic_ring:
        return False, "No aromatic ring found"

    # Look for amino group (-NH2 or -NH3+)
    amino_pattern = Chem.MolFromSmarts("[NH2,NH3]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"

    # Look for two-carbon chain connecting amino group and aromatic ring
    for atom_idx in amino_matches[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            return False, "Amino group is directly attached to the aromatic ring"
        neighbor_atoms = atom.GetNeighbors()
        for neighbor in neighbor_atoms:
            if neighbor.GetIsAromatic():
                continue
            path = Chem.FindAllPathsOfLengthN(mol, atom_idx, neighbor.GetIdx(), 3)
            if len(path) > 0:
                return True, "Monoamine structure detected"

    return False, "No two-carbon chain connecting amino group and aromatic ring"
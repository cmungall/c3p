"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""

from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    An N-hydroxy-alpha-amino-acid is an amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find alpha-amino acid backbone: N-C-C(=O)-O
    amino_acid_pattern = Chem.MolFromSmarts("N-C-C(=O)[O;H1]")
    aa_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not aa_matches:
        return False, "No alpha-amino acid backbone found"

    # For each match, check if nitrogen has at least one oxygen substituent
    for match in aa_matches:
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)

        # Get neighbors of nitrogen
        n_neighbors = nitrogen_atom.GetNeighbors()

        # Check for oxygen substituents on nitrogen
        has_n_oxygen = False
        for neighbor in n_neighbors:
            if neighbor.GetAtomicNum() == 8:
                has_n_oxygen = True
                break
            elif neighbor.GetAtomicNum() == 7:
                # Check for imino or oxime group (N=NOH)
                n2_neighbors = neighbor.GetNeighbors()
                for n2_neighbor in n2_neighbors:
                    if n2_neighbor.GetIdx() != nitrogen_idx and n2_neighbor.GetAtomicNum() == 8:
                        has_n_oxygen = True
                        break

        if has_n_oxygen:
            return True, "Contains N-hydroxy-alpha-amino-acid structure"

        # Check for double-bonded oxygen (N=O)
        for bond in nitrogen_atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(nitrogen_atom)
                if other_atom.GetAtomicNum() == 8:
                    has_n_oxygen = True
                    break

        if has_n_oxygen:
            return True, "Contains N-hydroxy-alpha-amino-acid structure"

    return False, "Does not contain N-hydroxy-alpha-amino-acid structure"
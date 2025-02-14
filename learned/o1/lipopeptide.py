"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with an attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bonds (amide bonds between amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C;!R]")  # Amide bond to non-ring carbon
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) < 3:
        return False, f"Insufficient peptide bonds found ({len(peptide_bonds)} found, need at least 3)"

    # Identify amino acid residues (N-C-C(=O) pattern)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H][CX3](=O)")  # N-C-C(=O)
    amino_acids = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acids) < 3:
        return False, f"Insufficient amino acid residues found ({len(amino_acids)} found, need at least 3)"

    # Function to find the longest aliphatic carbon chain
    def get_longest_aliphatic_chain(mol):
        max_length = 0
        chains = Chem.GetSymmSSSR(mol)
        for bond in mol.GetBonds():
            if bond.IsInRing():
                continue
            if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                path = Chem.FindAllPathsOfLengthN(mol, 8, useBonds=True, useHs=False)
                for p in path:
                    is_aliphatic = True
                    for idx in p:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetAtomicNum() != 6 or atom.IsInRing():
                            is_aliphatic = False
                            break
                    if is_aliphatic:
                        max_length = max(max_length, len(p))
        return max_length

    # Find the longest aliphatic chain
    longest_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            length = dfs_aliphatic_chain(atom, set(), mol)
            longest_chain_length = max(longest_chain_length, length)

    if longest_chain_length < 8:
        return False, f"No lipid chain found (longest aliphatic chain is {longest_chain_length} carbons, need at least 8)"

    # Verify that the lipid chain is connected to the peptide chain
    # Simplified check: molecule is connected
    if len(Chem.GetMolFrags(mol)) > 1:
        return False, "Molecule is disconnected"

    return True, "Contains both peptide bonds and lipid chain indicative of a lipopeptide"

def dfs_aliphatic_chain(atom, visited, mol):
    """
    Depth-first search to find the length of an aliphatic carbon chain.
    """
    length = 1
    visited.add(atom.GetIdx())
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing() and neighbor.GetIdx() not in visited:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() == Chem.BondType.SINGLE:
                length = max(length, 1 + dfs_aliphatic_chain(neighbor, visited.copy(), mol))
    return length
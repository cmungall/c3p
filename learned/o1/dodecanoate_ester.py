"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:36000 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the acyl chain (acid component) is derived from lauric acid (dodecanoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups: pattern [C](=O)O[C]
    ester_pattern = Chem.MolFromSmarts("[C:1](=O)[O:2][C,N]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check for lauric acid acyl chain
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[1]

        # Get the acyl chain length
        acyl_chain_length = get_acyl_chain_length(mol, carbonyl_c_idx, ester_o_idx)

        # Debugging output
        # print(f"Ester group at atoms {match}: acyl chain length = {acyl_chain_length}")

        if acyl_chain_length == 11:
            return True, "Contains dodecanoate ester group"

    return False, "No dodecanoate ester groups with lauric acid found"

def get_acyl_chain_length(mol, carbonyl_c_idx, ester_o_idx):
    """
    Traverses the acyl chain starting from the carbonyl carbon and counts the number of carbon atoms.

    Args:
        mol (RDKit Mol): The molecule
        carbonyl_c_idx (int): Index of the carbonyl carbon atom
        ester_o_idx (int): Index of the ester oxygen atom

    Returns:
        int: Number of carbon atoms in the acyl chain (excluding the carbonyl carbon)
    """
    visited = set()
    carbon_count = 0
    stack = []

    # Start from the neighbors of the carbonyl carbon, excluding the ester oxygen
    carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
    for neighbor in carbonyl_c_atom.GetNeighbors():
        n_idx = neighbor.GetIdx()
        if n_idx != ester_o_idx:
            stack.append(n_idx)

    while stack:
        idx = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)

        # Proceed only if the atom is a carbon atom
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(idx, neighbor_idx)
                # Traverse only single bonds to avoid double/triple bonds
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor_idx not in visited:
                    stack.append(neighbor_idx)
    return carbon_count
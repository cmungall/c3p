"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:35980 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the condensation of tetradecanoic acid (myristic acid) with an alcohol or phenol.
    It contains a linear, unbranched acyl chain of 14 carbons (including the carbonyl carbon) connected via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if the acyl chain is tetradecanoyl
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[1]     # Ester oxygen atom index

        # Traverse the acyl chain starting from the carbonyl carbon
        acyl_chain_atoms = set()
        visited = set([ester_o_idx])  # Exclude ester oxygen to avoid traversing into alcohol part

        def traverse_acyl_chain(atom_idx, prev_idx=None):
            if atom_idx in visited:
                return False  # Cycle detected
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                # Non-carbon atom, stop traversal
                return True
            acyl_chain_atoms.add(atom_idx)
            neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != prev_idx and n.GetIdx() != ester_o_idx]
            carbon_neighbors = [idx for idx in neighbors if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
            # Check for branching
            if len(carbon_neighbors) > 1:
                return False  # Branching detected
            for neighbor_idx in carbon_neighbors:
                if not traverse_acyl_chain(neighbor_idx, atom_idx):
                    return False
            return True  # Continue traversal

        if not traverse_acyl_chain(carbonyl_c_idx):
            continue  # Skip this ester group due to branching or cycle

        # Count the number of carbons in acyl chain (including carbonyl carbon)
        acyl_chain_length = len(acyl_chain_atoms)
        if acyl_chain_length == 14:
            return True, "Contains tetradecanoate ester group"

    return False, "No tetradecanoate ester groups found"
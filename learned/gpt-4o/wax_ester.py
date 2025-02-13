"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of a carboxy group of
    a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Get the atom indices of all ester linkages found
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Function to count contiguous carbon atoms using DFS
    def count_carbon_chain(start_atom_idx):
        visited = set()
        stack = [start_atom_idx]
        carbon_count = 0

        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # If carbon
                carbon_count += 1
                for neighbor in atom.GetNeighbors():
                    stack.append(neighbor.GetIdx())
        
        return carbon_count

    # Verify that each ester group connects two sufficiently long carbon chains
    for match in ester_matches:
        atom_idx_c1 = match[0]  # Carbon of the carbonyl group
        atom_idx_o = match[2]   # Oxygen in -C(=O)O-

        # Get neighbors excluding the ones directly involved in the ester linkage
        c1_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(atom_idx_c1).GetNeighbors() if n.GetIdx() != match[1]]
        o_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(atom_idx_o).GetNeighbors() if n.GetIdx() != match[1]]

        # Check if neighbors exist
        if not c1_neighbors or not o_neighbors:
            return False, "Missing chains connected to the ester linkage"

        c1_chain_length = max(count_carbon_chain(neigh_idx) for neigh_idx in c1_neighbors)
        o_chain_length = max(count_carbon_chain(neigh_idx) for neigh_idx in o_neighbors)

        # Ensure both chains are sufficiently long (e.g., >=8 carbons)
        if c1_chain_length < 8 or o_chain_length < 8:
            return False, f"Chains attached to ester group are too short (lengths found: {c1_chain_length} and {o_chain_length})"

    return True, "Contains an ester linkage with sufficiently long carbon chains typical of wax esters"
"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:24038 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Minimum acyl chain length for fatty acids (adjusted to 4)
    MIN_CHAIN_LENGTH = 4

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define amide functional group pattern (monocarboxylic acid amide)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")

    # Find all amide groups in the molecule
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amides = len(amide_matches)

    if num_amides == 0:
        return False, "No amide functional group found"

    if num_amides > 1:
        return False, f"Found {num_amides} amide groups; molecule is not a monocarboxylic acid amide"

    # Get the indices of the amide group
    amide_match = amide_matches[0]
    carbonyl_c_idx = amide_match[0]
    o_idx = amide_match[1]
    n_idx = amide_match[2]

    # Exclude carbonyl oxygen and amide nitrogen from traversal
    exclude_atoms = set([o_idx, n_idx])

    # Function to recursively traverse the acyl chain
    def traverse_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom_idx in visited:
            return 0
        if atom.IsInRing():
            return 0  # Exclude rings
        visited.add(atom_idx)
        if atom.GetAtomicNum() != 6:
            return 0  # Only follow carbons
        if atom.GetDegree() > 3:
            return 0  # Exclude branching
        length = 1  # Count current carbon
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in exclude_atoms]
        for nbr_idx in neighbors:
            if nbr_idx not in visited:
                length += traverse_chain(nbr_idx, visited.copy())
                break  # Only follow linear path (avoid branching)
        return length

    # Start traversal from the carbonyl carbon
    visited_atoms = set()
    acyl_chain_length = traverse_chain(carbonyl_c_idx, visited_atoms) - 1  # Exclude carbonyl carbon

    if acyl_chain_length >= MIN_CHAIN_LENGTH:
        return True, f"Molecule is a fatty amide with acyl chain length {acyl_chain_length}"
    else:
        return False, f"Acyl chain too short ({acyl_chain_length} carbons); not a fatty amide"

# Example usage:
# smiles = "CCCCCCCC(=O)NCCO"  # N-decanoylglycine
# result, reason = is_fatty_amide(smiles)
# print(result, reason)
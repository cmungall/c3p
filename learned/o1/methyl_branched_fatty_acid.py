"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclic structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic"

    # Identify carboxylic acid group using a specific SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 1"
    carboxyl_c_idx = matches[0][0]  # Index of the carboxyl carbon

    # Check for heteroatoms other than O in carboxyl group
    allowed_heteroatoms = set(matches[0])
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        idx = atom.GetIdx()
        if atomic_num == 6 or atomic_num == 1:
            continue  # Carbon or Hydrogen is allowed
        elif atomic_num == 8:
            if idx not in allowed_heteroatoms:
                return False, f"Extra oxygen atom found at index {idx + 1}"
        else:
            return False, f"Heteroatom {atom.GetSymbol()} found at index {idx + 1}"

    # Find the longest carbon chain starting from carboxyl carbon
    main_chain = find_longest_chain(mol, carboxyl_c_idx, set())

    # Check for branches (neighboring carbons not in main chain)
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() != 6:
                continue  # Skip non-carbon atoms
            if neighbor_idx in main_chain:
                continue  # Part of main chain
            # Found a branch, check if it's a methyl group
            branch_length = count_branch_length(mol, neighbor_idx, set([atom_idx]))
            if branch_length > 1:
                return False, f"Branch at atom {atom_idx + 1} is larger than a methyl group"

    return True, "Molecule is a methyl-branched fatty acid"

def find_longest_chain(mol, current_idx, visited):
    """
    Finds the longest carbon chain starting from the given atom index.

    Args:
        mol (Chem.Mol): RDKit molecule object
        current_idx (int): Current atom index
        visited (set): Set of visited atom indices

    Returns:
        list: List of atom indices representing the longest chain
    """
    visited.add(current_idx)
    atom = mol.GetAtomWithIdx(current_idx)
    max_chain = [current_idx]

    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx in visited:
            continue
        if neighbor.GetAtomicNum() != 6:
            continue  # Only traverse carbon atoms
        new_visited = visited.copy()
        sub_chain = find_longest_chain(mol, neighbor_idx, new_visited)
        if len(sub_chain) + 1 > len(max_chain):
            max_chain = [current_idx] + sub_chain

    return max_chain

def count_branch_length(mol, atom_idx, visited):
    """
    Counts the number of carbon atoms in a branch recursively.

    Args:
        mol (Chem.Mol): RDKit molecule object
        atom_idx (int): Starting atom index of the branch
        visited (set): Set of visited atom indices

    Returns:
        int: Number of carbon atoms in the branch
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 6:
        return 0
    visited.add(atom_idx)
    length = 1  # Count current carbon
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx in visited:
            continue
        if neighbor.GetAtomicNum() != 6:
            continue
        length += count_branch_length(mol, neighbor_idx, visited)
    return length
"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is a fatty acid in which the hydrocarbon chain has one or more alkyl substituents.
    The fatty acyl chain is usually saturated and the substituent is often a methyl group; however, unsaturated BCFAs 
    and branches other than methyl are also possible.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """

    from rdkit.Chem import rdmolops

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Assume the first carboxyl group
    carboxyl_carbon_idx = carboxylic_acid_matches[0][0]

    # Find all terminal carbons (carbons with only one carbon neighbor), excluding the carboxyl carbon
    terminal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() != carboxyl_carbon_idx:
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(neighbor_carbons) == 1:
                terminal_carbons.append(atom.GetIdx())

    if not terminal_carbons:
        return False, "No terminal carbons found in the hydrocarbon chain"

    # Find the longest carbon chain starting from the carboxyl carbon
    longest_chain = []
    max_length = 0
    for terminal_idx in terminal_carbons:
        # Find the shortest path between the terminal carbon and carboxyl carbon
        path = rdmolops.GetShortestPath(mol, carboxyl_carbon_idx, terminal_idx)
        # Verify that the path consists of carbons connected via single bonds
        valid_path = True
        for i in range(len(path) - 1):
            bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
            atom1 = mol.GetAtomWithIdx(path[i])
            atom2 = mol.GetAtomWithIdx(path[i+1])
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                valid_path = False
                break
            if atom1.GetAtomicNum() != 6 or atom2.GetAtomicNum() != 6:
                valid_path = False
                break
            if bond.IsInRing():
                valid_path = False
                break
        if valid_path and len(path) > max_length:
            max_length = len(path)
            longest_chain = path

    if not longest_chain:
        return False, "No suitable hydrocarbon chain found"

    # Check for branches
    branches = 0
    main_chain_set = set(longest_chain)
    for idx in longest_chain:
        atom = mol.GetAtomWithIdx(idx)
        neighbor_carbons = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        side_chain_atoms = [nbr_idx for nbr_idx in neighbor_carbons if nbr_idx not in main_chain_set]
        branches += len(side_chain_atoms)

    chain_length = len(longest_chain)
    if chain_length < 4:
        return False, f"Main chain length is {chain_length}, too short to be a fatty acid"

    if branches > 0:
        return True, f"Contains carboxyl group with hydrocarbon chain of length {chain_length} and {branches} branch(es)"
    else:
        return False, "No branches detected in the hydrocarbon chain"
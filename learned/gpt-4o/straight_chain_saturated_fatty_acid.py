"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if len(carboxylic_matches) != 1:
        return False, "No carboxylic acid group found or multiple groups present"
    
    # Check each atom in the carboxylic acid pattern to ensure it's at the end of a chain
    carbon, oxygen1, oxygen2 = carboxylic_matches[0]
    carbon_atom = mol.GetAtomWithIdx(carbon)

    # Get neighbors excluding carboxylic oxygens
    neighbors = [n.GetIdx() for n in carbon_atom.GetNeighbors() if n.GetIdx() not in {oxygen1, oxygen2}]
    if len(neighbors) != 1:
        return False, "Carboxylic acid group not at the end of a single continuous carbon chain"

    # Ensure no double or triple bonds in carbon chain, it should be saturated
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Presence of non-single C-C bonds indicates unsaturation"

    # Traverse the chain to ensure linear continuity and check for any branching
    visited = set()
    stack = [neighbors[0]]
    is_linear_chain = True
    hydroxy_count = 0

    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)

        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() == 8 and len([n for n in atom.GetNeighbors()]) == 1:
            hydroxy_count += 1
            if hydroxy_count > 1:
                return False, "More than one hydroxy group disrupting a straight chain"
        
        xposible_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(xposible_neighbors) > 2:
            is_linear_chain = False
            break
        
        stack.extend(xposible_neighbors)

    if not is_linear_chain:
        return False, "Branches or side chains detected; not a straight chain"

    return True, "Molecule is a straight-chain saturated fatty acid"
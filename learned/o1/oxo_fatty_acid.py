"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Assume only one carboxylic acid group for fatty acids
    carboxylic_acid_atom = None
    if len(carboxylic_acid_matches) == 1:
        # The first atom in the match is the carbonyl carbon
        carboxylic_acid_atom = carboxylic_acid_matches[0][0]
    else:
        return False, "Multiple carboxylic acid groups found"

    # Traverse the longest aliphatic chain from the carboxylic acid carbon
    visited = set()
    chain_atoms = []

    def traverse(atom_idx, chain_length, path):
        atom = mol.GetAtomWithIdx(atom_idx)
        visited.add(atom_idx)
        path.append(atom_idx)

        # Exclude atoms that are not carbon or are part of a ring
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            return 0, []

        max_length = chain_length
        max_path = list(path)

        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        for nbr_idx in neighbors:
            if nbr_idx not in visited:
                # Avoid going back to carboxylic acid oxygen atoms
                if nbr_idx in (carboxylic_acid_matches[0][1], carboxylic_acid_matches[0][2]):
                    continue
                length, nbr_path = traverse(nbr_idx, chain_length + 1, list(path))
                if length > max_length:
                    max_length = length
                    max_path = nbr_path

        return max_length, max_path

    # Start traversal from the carbon next to the carboxylic acid carbon
    max_chain_length = 0
    max_chain_path = []

    carboxylic_carbon = mol.GetAtomWithIdx(carboxylic_acid_atom)
    neighbors = [nbr.GetIdx() for nbr in carboxylic_carbon.GetNeighbors()]
    for nbr_idx in neighbors:
        # Skip the oxygens of the carboxylic acid
        if nbr_idx in (carboxylic_acid_matches[0][1], carboxylic_acid_matches[0][2]):
            continue
        visited.clear()
        length, path = traverse(nbr_idx, 1, [])
        if length > max_chain_length:
            max_chain_length = length
            max_chain_path = path

    # Check minimum chain length (e.g., at least 8 carbons)
    if max_chain_length < 8:
        return False, f"Aliphatic chain length is {max_chain_length}, which is too short for a fatty acid"

    # Exclude the carboxylic acid atoms from consideration in oxo group detection
    carboxylic_acid_atoms = set(carboxylic_acid_matches[0])

    # Check for additional aldehyde groups (-CHO)
    aldehyde = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = []
    for match in mol.GetSubstructMatches(aldehyde):
        if set(match).isdisjoint(carboxylic_acid_atoms):
            aldehyde_matches.append(match)

    # Check for ketone groups (>C=O)
    ketone = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = []
    for match in mol.GetSubstructMatches(ketone):
        if set(match).isdisjoint(carboxylic_acid_atoms):
            ketone_matches.append(match)

    total_oxo_groups = len(aldehyde_matches) + len(ketone_matches)

    if total_oxo_groups == 0:
        return False, "No additional aldehydic or ketonic group found"

    # Ensure that the oxo group is on the aliphatic chain
    oxo_on_chain = False
    chain_atom_set = set(max_chain_path)
    for match in aldehyde_matches + ketone_matches:
        if any(idx in chain_atom_set for idx in match):
            oxo_on_chain = True
            break

    if not oxo_on_chain:
        return False, "Oxo group is not on the aliphatic chain"

    return True, "Valid oxo fatty acid: long aliphatic chain with carboxylic acid and additional oxo group"
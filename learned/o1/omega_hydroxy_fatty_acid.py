"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid is a straight-chain fatty acid with a carboxyl group
    at position 1 and a hydroxyl group at the omega position (opposite end).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for rings - must be a linear molecule
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a straight-chain fatty acid"

    # Check for branches - all non-terminal carbons should have exactly two connections
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            degree = atom.GetDegree()
            idx = atom.GetIdx()
            if degree > 2:
                # Exclude terminal carbons (must have degree <= 3 in case of double bonds)
                neighbors = atom.GetNeighbors()
                if any(neigh.GetAtomicNum() != 1 and neigh.GetDegree() == 1 for neigh in neighbors):
                    continue  # Skip terminal carbons
                else:
                    return False, "Molecule is branched, not a straight-chain fatty acid"

    # Look for carboxyl group at one end
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Look for hydroxyl group at the omega end
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No terminal hydroxyl group found"

    # Confirm that hydroxyl group is at the end opposite the carboxyl group
    # Find the shortest path between the carboxyl carbon and the hydroxyl carbon
    min_path_length = None
    for carboxyl_match in carboxyl_matches:
        carboxyl_c_idx = carboxyl_match[0]
        for hydroxyl_match in hydroxyl_matches:
            hydroxyl_c_idx = hydroxyl_match[0]
            path = Chem.rdmolops.GetShortestPath(mol, carboxyl_c_idx, hydroxyl_c_idx)
            path_length = len(path) - 1  # Number of bonds between them
            if min_path_length is None or path_length < min_path_length:
                min_path_length = path_length
    if min_path_length is None:
        return False, "Cannot determine path between functional groups"

    # Check that the molecule is linear from carboxyl to hydroxyl group
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if min_path_length != num_carbons - 1:
        return False, "Hydroxyl group is not at the omega position"

    return True, "Molecule is an omega-hydroxy fatty acid"

__metadata__ = {
    'chemical_class': {
        'name': 'omega-hydroxy fatty acid',
        'definition': 'Any member of the class of naturally-occurring straight-chain fatty acids n carbon atoms long with a carboxyl group at position 1 and a hydroxyl at position n (omega).',
    }
}
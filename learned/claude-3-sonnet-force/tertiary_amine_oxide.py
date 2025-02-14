"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: CHEBI:36738 tertiary amine oxide

An N-oxide where there are three organic groups bonded to the nitrogen atom.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find N-oxide atoms
    n_oxides = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0)
    if n_oxides != 1:
        return False, f"Found {n_oxides} N-oxide groups, need exactly 1"

    # Find nitrogen atoms bonded to the N-oxide
    n_oxide_atom = next(atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0)
    n_atom = next(neighbor for neighbor in n_oxide_atom.GetNeighbors())

    # Check if nitrogen has 3 organic groups attached
    organic_neighbors = [neighbor for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() != 8]
    if len(organic_neighbors) != 3:
        return False, f"Found {len(organic_neighbors)} organic groups attached to N-oxide, need exactly 3"

    # Check if organic groups are not just hydrogen
    is_tertiary = all(atom.GetDegree() > 1 for atom in organic_neighbors)
    if not is_tertiary:
        return False, "One or more organic groups are just hydrogen"

    return True, "Contains an N-oxide with three organic groups attached"
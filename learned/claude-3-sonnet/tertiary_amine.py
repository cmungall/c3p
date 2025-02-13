"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:33899 tertiary amine
A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom with three substituents, all of which are carbon-containing groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if any nitrogen has exactly 3 carbon-containing substituents
    for nitrogen in nitrogen_atoms:
        carbon_substituents = [neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetAtomicNum() == 6]
        if len(carbon_substituents) == 3:
            return True, "Contains a nitrogen atom with 3 carbon-containing substituents"

    # No tertiary amine found
    return False, "No tertiary amine group found"
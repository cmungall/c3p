"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: CHEBI:58444 diketone
A compound that contains two ketone functionalities.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains exactly two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of ketone groups, excluding enol forms
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6](!@[OX2H1])")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    n_ketones = len(ketone_matches)

    if n_ketones == 2:
        return True, "Molecule contains exactly 2 ketone groups"
    else:
        return False, f"Found {n_ketones} ketone group(s), not a diketone"
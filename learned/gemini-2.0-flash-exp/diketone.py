"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a molecule containing exactly two ketone (C=O) functionalities.

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

    # Define ketone SMARTS pattern
    # [CX3](=[OX1])([CX4])([CX4])  is a carbon double bonded to oxygen and bonded to two other carbons
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])([CX4])([CX4])")

    # Find all matches of the ketone pattern
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Check if exactly 2 ketone groups are present
    if len(ketone_matches) == 2:
        return True, "Molecule contains exactly two ketone groups."
    else:
        return False, f"Molecule contains {len(ketone_matches)} ketone group(s), requires exactly two."
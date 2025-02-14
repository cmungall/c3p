"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is defined as any molecule containing exactly two carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OH1]")
    if carboxylic_acid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check number of matches
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_matches = len(matches)

    if num_matches == 2:
        return True, "Molecule contains exactly two carboxylic acid groups."
    elif num_matches < 2:
         return False, f"Molecule contains {num_matches} carboxylic acid groups, less than 2."
    else:
        return False, f"Molecule contains {num_matches} carboxylic acid groups, more than 2."
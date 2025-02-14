"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid, having a cyclohexane ring with multiple hydroxyl groups and a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the core quinic acid structure: a cyclohexane with 4 hydroxyls and 1 COOH
    quinic_acid_pattern = Chem.MolFromSmarts("[CX4]([OX2H1])([OX2H1])([OX2H1])[CX4]([CX4]([OX2H1])([OX2H1])[CX4])C(=O)[OX1H0]")
    if mol.HasSubstructMatch(quinic_acid_pattern):
        return True, "Matches the quinic acid pattern (cyclohexane ring, multiple hydroxyls, and a carboxyl group)"

    # Define a simpler pattern: a cyclohexane with at least 3 hydroxyl groups and 1 COOH group
    quinic_acid_pattern_simplified = Chem.MolFromSmarts("[CX4]1[CX4]([OX2H1])[CX4]([OX2H1])[CX4]([OX2H1])[CX4]([CX4]1)C(=O)[OX1H0]")
    if mol.HasSubstructMatch(quinic_acid_pattern_simplified):
        return True, "Matches simplified quinic acid pattern (cyclohexane ring with at least 3 hydroxyls and COOH)"

    return False, "Does not match the quinic acid pattern"
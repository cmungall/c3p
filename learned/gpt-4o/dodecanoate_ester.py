"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester should have lauric acid (12-carbon saturated chain) as the carboxylic acid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined lauric acid ester SMARTS pattern
    # This pattern is meant to specifically capture an ester bond with lauric acid
    lauric_acid_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O[CX4]")
    if lauric_acid_ester_pattern is None:
        return None, None  # If the pattern cannot be compiled

    # Check for lauric acid ester substructure
    if mol.HasSubstructMatch(lauric_acid_ester_pattern):
        return True, "Contains lauric acid ester component"

    return False, "Does not contain lauric acid ester component"
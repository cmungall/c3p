"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the esterification of tetradecanoic acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define tetradecanoate ester pattern (14-carbons with ester linkage)
    # [*] is a wildcard for any attachment, handles various alcohols
    tetradecanoate_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O[*]")
    
    # Check if the molecule matches the tetradecanoate ester pattern
    if mol.HasSubstructMatch(tetradecanoate_ester_pattern):
        return True, "Contains tetradecanoate ester group"
    else:
        return False, "No tetradecanoate ester group found"
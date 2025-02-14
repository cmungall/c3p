"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as any ester of carbamic acid or its N-substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a carbamate ester
    carbamate_ester_pattern = Chem.MolFromSmarts("[O][C](=[O])[N]")
    
    # Check if the molecule has the carbamate ester substructure
    if mol.HasSubstructMatch(carbamate_ester_pattern):
        return True, "Contains carbamate ester functional group"
    else:
        return False, "Does not contain carbamate ester functional group"

# Note: To use this function, pass a SMILES string as an argument, and it will return the classification.
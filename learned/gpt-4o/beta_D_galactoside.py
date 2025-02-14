"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside has a beta-configured D-galactose sugar moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-galactoside
    beta_d_galactoside_smarts = "[C@@H]1([C@@H]([C@@H](O)[C@H]([C@H](O1)CO)O)O)"

    # Create a molecule from the SMARTS pattern
    galactoside_pattern = Chem.MolFromSmarts(beta_d_galactoside_smarts)
    
    # Check if the molecule matches the beta-D-galactoside pattern
    if mol.HasSubstructMatch(galactoside_pattern):
        return True, "Contains beta-D-galactoside moiety"
    else:
        return False, "Beta-D-galactoside moiety not found"
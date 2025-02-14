"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is defined as a pyridoindole containing a beta-carboline skeleton
    and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-carboline core pattern (pyridoindole)
    # This pattern includes the core tricyclic structure of beta-carbolines 
    # assuming typical connectivity of pyridine and indole portions
    beta_carboline_pattern = Chem.MolFromSmarts("c1ccc2[nH]c3cccc(c3n2)c1")

    # Check if the molecule matches the beta-carboline pattern
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline skeleton"

    return False, "Does not contain beta-carboline skeleton"
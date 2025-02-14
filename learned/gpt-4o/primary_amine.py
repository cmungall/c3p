"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by one NH2 group attached to a carbon atom, not participating in peptide bonds (amides).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for primary amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][CX4]")
    
    # Check for matches of the primary amine pattern
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Primary amine group (NH2) attached to carbon detected"

    return False, "No primary amine group (NH2) attached to carbon detected"
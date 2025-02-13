"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is characterized by the presence of a five-membered lactone ring; various substituents and 
    structural variations such as double bonds or stereochemistry centers are considered.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a gamma-lactone, False otherwise.
        str: Explanation for the classification decision.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a comprehensive SMARTS pattern for gamma-lactone structure
    # Pattern is flexible to allow various substitutions and structural flexibilities
    gamma_lactone_patterns = [
        Chem.MolFromSmarts("O=C1OCC=C1"),   # Open double bond lactone
        Chem.MolFromSmarts("O=C1OCCC1"),    # Basic lactone
        Chem.MolFromSmarts("O=C1OC=C[CH2]1"), # Lactone with alternative ring saturation
        Chem.MolFromSmarts("O=C1OC(=O)C=C1"),  # Highly oxidized lactone
        Chem.MolFromSmarts("O=C1OC(C)=CC1")  # Substituted gamma lactone
        # More patterns can be added if other configurations are frequently found
    ]

    for pattern in gamma_lactone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a five-membered lactone ring characteristic of gamma-lactones"

    return False, "No five-membered lactone ring found"
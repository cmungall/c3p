"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a versatile gamma-lactone pattern for 2-furanone which allows for substituents
    gamma_lactone_patterns = [
        Chem.MolFromSmarts("O=C1OCC=C1"),  # Basic 2-furanone core
        Chem.MolFromSmarts("O=C1O/C=C/C1"),  # Tolerates other bond changes
        Chem.MolFromSmarts("O=C1[C;R1]=[C;R1]OC1"),  # Tolerate substitution variations
    ]
    
    # Check if the structure matches any of the gamma-lactone patterns
    for pattern in gamma_lactone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Structure matches gamma-lactone (2-furanone) backbone pattern"
    
    return False, "No gamma-lactone (2-furanone) pattern found"

# Example usage
smiles = "O=C1OCC=C1"  # Example of butenolide (but-2-en-4-olide)
result, reason = is_butenolide(smiles)
print(result, reason)
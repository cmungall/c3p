"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_6_aminopurine(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a 6-aminopurine moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define 6-aminopurine (adenine) pattern
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    
    # Check if molecule contains adenine pattern
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) moiety"
    else:
        return False, "Does not contain 6-aminopurine (adenine) moiety"

# Example usage
smiles1 = "Nc1ncnc2n(cnc12)C3OC(COP(=O)(O)OP(=O)(O)OCC4OC(n5cnc6c(N)ncnc56)C(O)C4O)C(O)C3O"  # Acetyl-CoA
smiles2 = "CCc1ccccc1"  # Ethylbenzene

print(is_6_aminopurine(smiles1))  # Output: (True, 'Contains 6-aminopurine (adenine) moiety')
print(is_6_aminopurine(smiles2))  # Output: (False, 'Does not contain 6-aminopurine (adenine) moiety')
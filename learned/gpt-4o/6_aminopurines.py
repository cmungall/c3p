"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) structure based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains the 6-aminopurine structure, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 6-aminopurine (adenine) SMARTS pattern
    # Adenine's SMILES representation is: "c1ncnc2ncnc12"
    # Adding the 6-amino group: "c1ncnc2nc(N)nc12"
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2nc(N)nc12")
    if adenine_pattern is None:
        return False, "Failed to create 6-aminopurine pattern"

    # Check if the molecule contains the adenine structure
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) structure"
    else:
        return False, "Does not contain 6-aminopurine (adenine) structure"

# Example usage for demonstration purposes
# smiles_example = "C1=NC2=C(N1)N=C(NC2=N)N"  # SMILES for adenine
# print(is_6_aminopurines(smiles_example))
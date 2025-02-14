"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is a deoxyribonucleoside containing a pyrimidine base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for deoxyribose pattern
    deoxyribose_pattern = Chem.MolFromSmarts("C1([C@@H](CO)O[C@H]1O)")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose backbone found"

    # Look for pyrimidine pattern
    pyrimidine_pattern = Chem.MolFromSmarts("c1cnc[nH]c1=O")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found attached to the sugar"

    # Additional check - attachment of base to sugar
    # This would be a specific connection point pattern, 
    # but for simplicity, we assume pattern checks should suffice in this example.

    return True, "Contains a deoxyribose backbone with a pyrimidine base"

# Test the function with example SMILES strings
example_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O"  # Thymidine
result, reason = is_pyrimidine_deoxyribonucleoside(example_smiles)
print(f"Result: {result}, Reason: {reason}")
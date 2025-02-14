"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine group is characterized by the presence of a -C(=NR)-NR2 group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a carboxamidine group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS pattern for carboxamidine: Ensure it captures the full characteristic group
    carboxamidine_pattern = Chem.MolFromSmarts("C(=N)N")

    # Check for the presence of the refined carboxamidine pattern
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine group (-C(=NR)NR2)"

    return False, "Does not contain carboxamidine group"

# Example usage:
# test_smiles = "NC(=N)c1ccccc1"  # benzamidine
# print(is_carboxamidine(test_smiles))
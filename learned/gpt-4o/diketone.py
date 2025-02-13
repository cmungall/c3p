"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone (C=O) functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a ketone group, considering potential ring structures
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")

    # Find ketone groups in the molecule, making sure realistic C=O occurrences are found
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Filter out duplicates in cyclic ketone environments by checking atom indices
    unique_ketone_bonds = set((min(m[0], m[1]), max(m[0], m[1])) for m in ketone_matches)

    # Check if there are exactly two ketone groups
    if len(unique_ketone_bonds) == 2:
        return True, "Contains exactly two ketone groups, hence a diketone."
    else:
        return False, f"Contains {len(unique_ketone_bonds)} ketone groups, not a diketone."

# Examples
# Test the function with a diketone example SMILES
result, reason = is_diketone("CCCCCC(=O)CC(=O)CCCCCC")  # Expected: True, diketone
print(result, reason)

# Test the function with a non-diketone example SMILES
result, reason = is_diketone("CCCCCCC=O")  # Expected: False, only one ketone
print(result, reason)
"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: CHEBI:36630 diketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains two ketone functionalities.

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
    
    # Count ketone groups
    ketone_pattern = Chem.MolFromSmarts("C(=O)[CH3]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    n_ketones = len(ketone_matches)
    
    if n_ketones == 2:
        return True, "Contains two ketone functionalities"
    elif n_ketones < 2:
        return False, f"Found only {n_ketones} ketone groups, need two"
    else:
        return False, "Found more than two ketone groups"

# Example usage
smiles = "CC(=O)CC(=O)C"  # 2,4-Pentanedione
is_diketone, reason = is_diketone(smiles)
print(f"Is diketone? {is_diketone} ({reason})")
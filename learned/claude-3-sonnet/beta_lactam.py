"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:25863 beta-lactam
A lactam in which the amide bond is contained within a four-membered ring, which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-lactam ring pattern:
    # 4-membered ring with nitrogen, carbonyl carbon, and 2 other atoms (usually carbons)
    beta_lactam_pattern = Chem.MolFromSmarts("[NR1]1[CR1][CR1][CR1]C(=O)1")
    beta_lactam_matches = mol.GetSubstructMatches(beta_lactam_pattern)
    
    if beta_lactam_matches:
        return True, "Contains a beta-lactam ring (4-membered ring with N, C=O, and 2 other atoms)"
    else:
        return False, "No beta-lactam ring found"

# Example usage:
print(is_beta_lactam("CC1(C(N2C(S1)C(C2=O)NC(=O)C(N)C3=CC=C(O)C=C3)C(=O)O)C"))  # True
print(is_beta_lactam("CCOC(=O)C1=CC=CC=C1"))  # False
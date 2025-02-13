"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring (a four-membered cyclic amide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for beta-lactam ring pattern (four-membered ring with a C=O and a nitrogen)
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)NC1")
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a beta-lactam ring"
    else:
        return False, "No beta-lactam ring found"

# Example usage:
# result, reason = is_beta_lactam_antibiotic('CC1=C(C(=NN1CC(=O)NC2C3N(C2=O)C(=C(CS3)CSC4=NN=C(S4)C)C(=O)O)C)Cl')
# print(result, reason)
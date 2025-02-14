"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:27731 beta-lactam antibiotic
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

    # Define SMARTS pattern for the beta-lactam ring (four-membered cyclic amide)
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N[C]1")

    # Check for substructure match
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a beta-lactam ring"
    else:
        return False, "Does not contain a beta-lactam ring"
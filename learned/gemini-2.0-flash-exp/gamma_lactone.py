"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a lactone with a five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the gamma-lactone ring
    gamma_lactone_pattern = Chem.MolFromSmarts('[C,c]1[C,c][C,c](=O)[O][C,c]1')

    if gamma_lactone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if the molecule has the gamma-lactone pattern
    matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    if not matches:
          return False, "No gamma-lactone ring pattern found."
    
    return True, "Molecule contains a gamma-lactone ring."
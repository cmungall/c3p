"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 2-furanone skeleton (gamma-lactone with double bond)
    butenolide_pattern = Chem.MolFromSmarts('O=C1C=CCO1')
    if butenolide_pattern is None:
        return False, "Error in SMARTS pattern for butenolide"

    # Check for the 2-furanone skeleton
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains 2-furanone skeleton"
    else:
        return False, "Does not contain 2-furanone skeleton"
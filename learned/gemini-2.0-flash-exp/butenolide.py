"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone, a five-membered ring containing an ester group
    with the carbonyl carbon adjacent to the oxygen in the ring.

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

    # Define the SMARTS pattern for the butenolide core (2-furanone)
    butenolide_pattern = Chem.MolFromSmarts("[OX2]1[CX3](=[OX1])[CX3]=[CX3]1")

    # Check if the molecule has the butenolide core structure
    if mol.HasSubstructMatch(butenolide_pattern):
         return True, "Has butenolide core structure (2-furanone)"
    else:
        return False, "Does not have butenolide core structure (2-furanone)"
"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide
    to the C=O bond of an aldehyde or ketone.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for cyanohydrin functional group
    cyanohydrin_pattern = Chem.MolFromSmarts("[C;X4]([OH])([C]#[N])")
    if cyanohydrin_pattern is None:
        return False, "Invalid SMARTS pattern for cyanohydrin"

    # Search for the cyanohydrin substructure
    matches = mol.GetSubstructMatches(cyanohydrin_pattern)
    if not matches:
        return False, "Cyanohydrin functional group not found"

    # If matches are found, return True with a reason
    return True, "Contains the cyanohydrin functional group (alpha-hydroxynitrile)"
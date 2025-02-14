"""
Classifies: CHEBI:23437 cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carbon with attached hydroxyl group (OH) and nitrile group (C#N)
    # Exclude patterns that might encompass ring structures
    cyanohydrin_pattern = Chem.MolFromSmarts("[C;!R]([OH])([C]#[N])")

    # Check if molecule matches the cyanohydrin pattern
    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Molecule contains the defining substructure of a cyanohydrin"
    else:
        return False, "No cyanohydrin substructure found in the molecule"
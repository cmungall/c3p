"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a
    saturated carbon atom which has either three hydrogen atoms attached to it or
    only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for primary alcohol carbons, matching a -CH2-OH group 
    # with the carbon being saturated, not bonded to any other non-hydrogen, non-oxygen
    # atom by a double or triple bond
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2;!$(C=[C,N,O,S]);!$(C#[C,N])][OH]")

    # Find matches for the pattern
    matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    
    # Check if any of the patterns are present
    if not matches:
        return False, "No primary alcohol group found"

    return True, "Primary alcohol group found"
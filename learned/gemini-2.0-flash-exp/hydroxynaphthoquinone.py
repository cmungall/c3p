"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone with at least one hydroxy substituent.
    The naphthoquinone must have carbonyl groups in 1-4 relative positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for naphthoquinone core and hydroxy
    # Using SMARTS that checks for a 1-4 naphthoquinone.
    # This pattern uses query atoms 'a' to represent an aromatic carbon
    naphthoquinone_pattern = Chem.MolFromSmarts("a12ccccc1C(=O)C=CC2=O")
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")

    # Check if core is present
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core found"

    # Check if at least one hydroxy is present
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group found"

    return True, "Contains a naphthoquinone core and at least one hydroxy substituent"
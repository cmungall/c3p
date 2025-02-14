"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for pyrrole ring
    # Matches an aromatic five-membered ring with one nitrogen (bearing a hydrogen) and four carbons
    pyrrole_smarts = '[nH]1ccccc1'  # Incorrect, actually includes a six-membered ring. We need to correct this.

    # Correct SMARTS pattern for pyrrole
    pyrrole_smarts = '[nH]1cccc1'  # Correct pattern for pyrrole ring
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)

    # Find all matches of the pyrrole substructure
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    pyrrole_count = len(pyrrole_matches)

    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole units"
    else:
        return False, f"Contains only {pyrrole_count} pyrrole unit(s), less than 2 required for polypyrrole"
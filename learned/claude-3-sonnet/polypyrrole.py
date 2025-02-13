"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: CHEBI:38575 polypyrrole
A compound composed of two or more pyrrole units.
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Create a pyrrole ring query
    pyrrole_query = rdqueries.IsNPyrroles(2)

    # Check if molecule contains at least two pyrrole rings
    for match in mol.GetSubstructMatches(pyrrole_query, maxMatches=2):
        if len(match) == 2:
            return True, "Contains at least two pyrrole units"

    return False, "Does not contain at least two pyrrole units"
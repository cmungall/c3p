"""
Classifies: CHEBI:38077 polypyrrole
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

    # Define the pyrrole substructure using SMARTS (relaxed version)
    pyrrole_smarts = "[n;H0,H1]1cc[c,n]c1"
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)

    # Get all matches of the pattern
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)

    # Count the number of unique pyrrole rings by counting matches that dont overlap by an edge
    # Initialize a set to track atoms that are in a pyrrole ring
    matched_atoms = set()
    num_pyrroles = 0

    # Loop through all matches.
    for match in pyrrole_matches:
        #Convert the tuple to a set of integers so we can check for overlaps
        match_set = set(match)
        # Check to see if there are any atoms shared between the current pyrrole, and all the other pyrroles seen so far.
        if len(matched_atoms.intersection(match_set)) == 0:
            num_pyrroles += 1
            matched_atoms.update(match_set)

    # Classify based on the count of pyrrole rings
    if num_pyrroles >= 2:
        return True, f"Contains {num_pyrroles} pyrrole units, therefore a polypyrrole."
    else:
        return False, f"Contains only {num_pyrroles} pyrrole unit(s), therefore not a polypyrrole."
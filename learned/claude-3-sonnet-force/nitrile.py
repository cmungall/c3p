"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:35698 nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is a compound having the structure RC#N.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C#N]")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group (-C#N) found"

    # Count number of nitrile groups
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    n_nitriles = len(nitrile_matches)

    # Nitriles should have one or more nitrile groups
    if n_nitriles == 0:
        return False, "No nitrile groups found"

    return True, f"Molecule contains {n_nitriles} nitrile group(s)"
"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is a compound having the structure R-C#N, where 'R' can be any substituent except hydrogen.
    The suffix nitrile denotes the triply bound #N atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define allowed elements (common in organic molecules)
    allowed_atomic_nums = set([1, 6, 7, 8, 9, 15, 16, 17, 35, 53])  # H, C, N, O, F, P, S, Cl, Br, I

    # Check for disallowed elements (metals, etc.)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Contains metal or disallowed atoms; likely not an organic nitrile"

    # Define the nitrile SMARTS pattern:
    # [!#1][C]#[N]: Nitrile carbon connected to any atom except hydrogen
    nitrile_pattern = Chem.MolFromSmarts("[!#1][C]#[N]")
    if nitrile_pattern is None:
        return False, "Invalid nitrile SMARTS pattern"

    # Search for nitrile substructure
    if mol.HasSubstructMatch(nitrile_pattern):
        return True, "Contains nitrile group (R-C#N) per definition"
    else:
        return False, "Does not contain nitrile group (R-C#N) per definition"
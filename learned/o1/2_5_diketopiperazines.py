"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazine
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is defined as any piperazinone that has a piperazine-2,5-dione skeleton,
    which is a six-membered ring with nitrogen atoms at positions 1 and 4, and carbonyl groups at
    positions 2 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the 2,5-diketopiperazine core
    diketopiperazine_pattern = Chem.MolFromSmarts('O=C1NCC(=O)NC1')  # piperazine-2,5-dione skeleton
    if diketopiperazine_pattern is None:
        return None, "Invalid SMARTS pattern for diketopiperazine"

    # Check if the molecule contains the 2,5-diketopiperazine core
    if mol.HasSubstructMatch(diketopiperazine_pattern):
        return True, "Molecule contains piperazine-2,5-dione skeleton"
    else:
        return False, "Molecule does not contain piperazine-2,5-dione skeleton"
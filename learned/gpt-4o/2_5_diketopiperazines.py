"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton, featuring a six-membered ring
    with alternating carbon and nitrogen atoms and ketone groups at positions 2 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a 2,5-diketopiperazine
    diketopiperazine_pattern = Chem.MolFromSmarts("C1(=O)NCC(=O)NC1")
    if mol.HasSubstructMatch(diketopiperazine_pattern):
        return True, "Contains a piperazine-2,5-dione skeleton"
    
    return False, "Does not contain a piperazine-2,5-dione skeleton"
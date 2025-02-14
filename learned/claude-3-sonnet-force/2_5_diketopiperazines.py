"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:34614 2,5-diketopiperazine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_5_diketopiperazine(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is a piperazine-2,5-dione skeleton.

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

    # Look for piperazine-2,5-dione core
    core_pattern = Chem.MolFromSmarts("O=C1NC(=O)CN1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No piperazine-2,5-dione core found"
    
    # Check for substituents on piperazine ring
    substituted_pattern = Chem.MolFromSmarts("C1NC(=O)CN(C)C1=O")
    if mol.HasSubstructMatch(substituted_pattern):
        return True, "Contains piperazine-2,5-dione core with substituents"
    else:
        return True, "Contains unsubstituted piperazine-2,5-dione core"

    # If none of the above conditions are met, return False
    return False, "Unknown reason"
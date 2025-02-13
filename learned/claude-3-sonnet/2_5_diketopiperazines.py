"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:50813 2,5-diketopiperazines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_5_diketopiperazine(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine contains a piperazine-2,5-dione skeleton.

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

    # Look for piperazine-2,5-dione skeleton pattern
    pattern = Chem.MolFromSmarts("O=C1NC(=O)CN=C1")
    if not mol.HasSubstructMatch(pattern):
        return False, "No piperazine-2,5-dione skeleton found"

    # Check for cyclic structure
    if not mol.GetRingInfo().IsAtomRingBond(0):
        return False, "Structure is not cyclic"

    # Count nitrogens and oxygens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count != 2:
        return False, "Must have exactly 2 nitrogen atoms"
    if o_count != 2:
        return False, "Must have exactly 2 oxygen atoms"

    return True, "Contains a piperazine-2,5-dione skeleton"
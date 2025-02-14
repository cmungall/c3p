"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:16581 mononitrophenol

A mononitrophenol is a phenol carrying a single nitro substituent at an unspecified position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phenol ring
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol ring found"
    
    # Check for single nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, expected 1"
    
    # Check for aromaticity
    if not mol.GetAromaticRings():
        return False, "Not aromatic"
    
    # Count carbons, hydrogens, nitrogens and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6 or c_count > 12:
        return False, "Wrong number of carbons for mononitrophenol"
    if h_count < 5 or h_count > 9:
        return False, "Wrong number of hydrogens for mononitrophenol"
    if n_count != 1:
        return False, "Expected exactly 1 nitrogen for mononitrophenol"
    if o_count != 2:
        return False, "Expected exactly 2 oxygens for mononitrophenol"
    
    return True, "Contains a phenol ring with a single nitro substituent"
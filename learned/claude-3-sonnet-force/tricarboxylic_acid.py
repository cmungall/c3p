"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:33498 tricarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as an oxoacid containing three carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxy group pattern (-C(=O)O)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Check if there are exactly 3 carboxy groups
    if len(carboxy_matches) != 3:
        return False, f"Found {len(carboxy_matches)} carboxy groups, need exactly 3"
    
    # Check for oxoacid (contains at least one C=O group)
    oxo_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo (C=O) group found, not an oxoacid"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 4:
        return False, "Too few carbons for tricarboxylic acid"
    if o_count < 6:
        return False, "Too few oxygens for tricarboxylic acid"
    
    return True, "Contains three carboxy groups and is an oxoacid"
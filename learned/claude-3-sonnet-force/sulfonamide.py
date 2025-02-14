"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:24062 sulfonamide
An amide of a sulfonic acid RS(=O)2NR'2, where R' can be H, alkyl, aryl, etc.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonamide groups
    sulfonamide_pattern1 = Chem.MolFromSmarts("S(=O)(=O)N")  # General sulfonamide pattern
    sulfonamide_pattern2 = Chem.MolFromSmarts("S(=O)(=O)NC")  # Sulfonamide with C attached to N
    sulfonamide_pattern3 = Chem.MolFromSmarts("S(=O)(=O)N([H])")  # Sulfonamide with H attached to N
    
    # Count sulfonamide groups (allow multiple)
    sulfonamide_matches = mol.GetSubstructMatches(sulfonamide_pattern1) \
                        + mol.GetSubstructMatches(sulfonamide_pattern2) \
                        + mol.GetSubstructMatches(sulfonamide_pattern3)
    
    if not sulfonamide_matches:
        return False, "No sulfonamide groups found"
    
    # Check for common sulfonamide substructures
    aromatic_sulfonamide_pattern = Chem.MolFromSmarts("c1ccccc1S(=O)(=O)N")
    aliphatic_sulfonamide_pattern = Chem.MolFromSmarts("CCS(=O)(=O)N")
    
    has_aromatic_sulfonamide = mol.HasSubstructMatch(aromatic_sulfonamide_pattern)
    has_aliphatic_sulfonamide = mol.HasSubstructMatch(aliphatic_sulfonamide_pattern)
    
    # Avoid false positives with sulfone groups
    sulfone_pattern = Chem.MolFromSmarts("S(=O)(=O)C")
    if mol.HasSubstructMatch(sulfone_pattern):
        return False, "Contains sulfone groups, not a sulfonamide"
    
    if has_aromatic_sulfonamide or has_aliphatic_sulfonamide:
        return True, "Contains one or more sulfonamide groups"
    else:
        return False, "Sulfonamide group not in expected environment"
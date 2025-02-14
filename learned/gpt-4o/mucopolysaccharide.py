"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely to be a mucopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has characteristics of a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Simplified check for uronic acid-like and glycosamine-like features
    # Uronic acids (often COOH groups), glycosamine (often NH2 groups)
    uronic_acid_pattern = Chem.MolFromSmarts("[C,c](=O)[O,N]")
    glycosamine_pattern = Chem.MolFromSmarts("[N;H2]")
    
    uronic_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    
    if len(uronic_matches) == 0 or len(glycosamine_matches) == 0:
        return False, "Missing uronic acid or glycosamine features"
    
    # Assume partial esterification with sulfate is present more heuristically,
    # as strict identification would require specific sulfate group context
    sulfate_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX2])")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    
    # Simplified logic: detect presence of all components, not ratios or polymer length
    if len(sulfate_matches) > 0:
        return True, "Characteristics of mucopolysaccharide present, including sulfate groups"
    else:
        return True, "Characteristics of mucopolysaccharide present, but no sulfate groups"

    return False, "Could not classify as a mucopolysaccharide"
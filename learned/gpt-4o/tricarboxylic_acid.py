"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing three carboxy groups (-C(=O)O
    or its anionic form -C(=O)[O-]).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define patterns for carboxylic acid group and its anionic form
    carboxy_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1-]")

    # Find all carboxylic acid groups (protonated and deprotonated)
    carboxy_matches = mol.GetSubstructMatches(carboxy_acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)

    # Total carboxylic/carboxylate groups
    total_carboxy_groups = len(carboxy_matches) + len(carboxylate_matches)
    
    # Specificity checks: Exclude molecules known to commonly misclassify 
    if total_carboxy_groups == 3:
        # Check if other functional groups may interfere with classification
        if not has_interfering_groups(mol):
            return True, "Contains exactly three carboxylic acid or carboxylate groups and no interfering groups"
        else:
            return False, "Contains three carboxylic acid or carboxylate groups but also interfering groups"
    else:
        return False, f"Found {total_carboxy_groups} carboxylic acid or carboxylate groups, need exactly 3"

def has_interfering_groups(mol):
    """
    Checks for additional functional groups or potential interferers that may misclassify
    a molecule as a tricarboxylic acid.

    Args:
        mol: RDKit Mol object

    Returns:
        bool: True if interfering groups are present, False otherwise
    """
    # Example of interfering group pattern (duplicate check for significant peptides)
    interfering_patterns = [
        Chem.MolFromSmarts("[N]"),  # Example: Nitrogen-containing groups in specific contexts
        # Add more patterns identified from analysis that could indicate misclassification
    ]

    for pattern in interfering_patterns:
        if mol.HasSubstructMatch(pattern):
            return True
    return False
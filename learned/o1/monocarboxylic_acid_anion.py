"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all carboxyl groups (protonated or deprotonated)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[O-,OH]")
    carboxyl_groups = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl_groups = len(carboxyl_groups)
    
    if num_carboxyl_groups != 1:
        return False, f"Found {num_carboxyl_groups} carboxyl groups, expected exactly 1"
    
    # Find deprotonated carboxylate groups
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_groups = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate_groups = len(carboxylate_groups)
    
    if num_carboxylate_groups != 1:
        return False, "Carboxyl group is not deprotonated"
    
    # If all conditions are met, classify as monocarboxylic acid anion
    return True, "Molecule is a monocarboxylic acid anion with exactly one deprotonated carboxyl group"
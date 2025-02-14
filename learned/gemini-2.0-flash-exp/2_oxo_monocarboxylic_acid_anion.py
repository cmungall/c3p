"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion has a C=O group at the 2-position relative to a carboxylic acid anion (C(=O)[O-]) group.
    It has exactly one carboxylic acid group

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    #Define 2-oxo monocarboxylic acid anion pattern
    oxo_acid_pattern = Chem.MolFromSmarts("C(=O)C(=O)[O-]")
    
    # Define carboyxlic acid pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O-]")

    # Check if molecule matches the 2-oxo monocarboxylic acid anion pattern
    if not mol.HasSubstructMatch(oxo_acid_pattern):
        return False, "Molecule does not contain the required 2-oxo monocarboxylic acid anion substructure"

    # Check if molecule has exactly one carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Molecule contains {len(carboxylic_acid_matches)} carboxylic acid groups, expected exactly 1"
    

    return True, "Molecule is a 2-oxo monocarboxylic acid anion"
"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-Oxo Monocarboxylic Acid Anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion
    based on its SMILES string.

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

    # Define SMARTS pattern for carboxylate group with adjacent oxo group
    oxo_monocarboxylic_pattern = Chem.MolFromSmarts("C(=O)C([O-])=O")

    # Check if the pattern matches within the molecule
    if mol.HasSubstructMatch(oxo_monocarboxylic_pattern):
        return True, "Contains an oxo group in the 2-position adjacent to a carboxylate group"
    
    return False, "Does not contain a 2-oxo group adjacent to a carboxylate group"
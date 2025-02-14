"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: ChEBI:51624 2-oxo monocarboxylic acid
Any monocarboxylic acid having a 2-oxo substituent.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for a 2-oxo group (-C(=O)-C-)
    oxo_pattern = Chem.MolFromSmarts("C(=O)C")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No 2-oxo group found"
    
    # Check if the 2-oxo group is attached to the carboxylic acid group
    carboxylic_acid_atoms = mol.GetSubstructMatch(carboxylic_acid_pattern)
    for oxo_match in oxo_matches:
        if oxo_match[1] in carboxylic_acid_atoms:
            return True, "Contains a 2-oxo monocarboxylic acid group"
    
    # If no 2-oxo group is attached to the carboxylic acid group, it's not a 2-oxo monocarboxylic acid
    return False, "The 2-oxo group is not attached to the carboxylic acid group"
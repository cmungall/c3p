"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide consists of four monosaccharide units linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to identify common monosaccharide units (furanose/pyranose rings with OH groups)
    # This pattern needs to capture variations of such units realistically
    monosaccharide_pattern = Chem.MolFromSmarts("O[C@H]1(CO)[C@@H]([C@H](O)[C@@H](O)O1) | [C@H]1(O[C@H]([C@H]([C@@H](O)C1O)O)CO)[C@H](O)O")
    
    if not monosaccharide_pattern:
        return (None, None)
    
    # Find all matches of the pattern in the molecule
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    # Count the number of unique monosaccharide units (this might need adjustment depending on exact determination)
    num_monosaccharide_units = len(monosaccharide_matches)
    
    if num_monosaccharide_units == 4:
        return True, "Contains four monosaccharide units linked together"
    else:
        return False, f"Found {num_monosaccharide_units} monosaccharide units, need exactly 4"

# Example usage:
# print(is_tetrasaccharide("OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)COC[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"))  # Modify with example SMILES
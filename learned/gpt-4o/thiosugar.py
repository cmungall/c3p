"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has one or more oxygens or hydroxy groups replaced by sulfur or sulfur-containing groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sugar framework (e.g., pyran, furan rings with hydroxyl groups)
    # In SMARTS representation, this would be a cyclic ether with hydroxyl (O) groups
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H]([C@H](O)[C@H](O)[C@H]1O)O)CO")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar backbone found"

    # Check for sulfur substitution (-S- or -SH)
    sulfur_pattern = Chem.MolFromSmarts("[#16]")  # Sulfur atom
    if not mol.HasSubstructMatch(sulfur_pattern):
        return False, "No sulfur substitution found"
    
    # Check for cases where sulfur atom is attached in place of a hydroxy group
    thio_pattern = Chem.MolFromSmarts("[C@@H](-[S])")
    if not mol.HasSubstructMatch(thio_pattern):
        return False, "No thio substitution in sugar structure"

    return True, "Contains sugar backbone with sulfur substitution"
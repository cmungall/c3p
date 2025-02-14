"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:49144 carbapenems
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are a class of beta-lactam antibiotics with a carbapenem skeleton
    that is variously substituted at positions 3, 4, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbapenem backbone pattern
    backbone_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2CC(=C(N1C(=O)C2)[C@@H](O)C)S")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No carbapenem backbone found"
    
    # Look for beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2C(=O)N1C2")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    # Check for substituents at positions 3, 4, and 6
    sub_pattern_3 = Chem.MolFromSmarts("[C@@H](C)[C@@H](O)[C@@H]1C(=O)N2C(=C1[C@H](C)[C@H]2S)")
    sub_pattern_4 = Chem.MolFromSmarts("[C@@H](C)[C@@H](O)[C@@H]1C(=O)N2C(=C1[C@H](C)[C@@H]2S)")
    sub_pattern_6 = Chem.MolFromSmarts("[C@@H](C)[C@@H](O)[C@@H]1C(=O)N2C(=C1[C@H](C)[C@H]2S)")
    
    if not (mol.HasSubstructMatch(sub_pattern_3) or mol.HasSubstructMatch(sub_pattern_4) or mol.HasSubstructMatch(sub_pattern_6)):
        return False, "No substitution at positions 3, 4, or 6 on carbapenem backbone"
    
    # Check for common carbapenem functional groups
    sulfur_pattern = Chem.MolFromSmarts("S")
    amine_pattern = Chem.MolFromSmarts("N")
    hydroxy_pattern = Chem.MolFromSmarts("O[H]")
    
    if not (mol.HasSubstructMatch(sulfur_pattern) and mol.HasSubstructMatch(amine_pattern) and mol.HasSubstructMatch(hydroxy_pattern)):
        return False, "Missing common functional groups (sulfur, amine, hydroxyl)"
    
    return True, "Contains carbapenem skeleton with substitutions at positions 3, 4, and/or 6"
"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMILES patterns for monosaccharide rings and glycosidic linkages
    monosaccharide_pattern = Chem.MolFromSmarts("[C@H]1(O)O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1")
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@H]")
    
    # Check for at least two monosaccharide units
    if len(mol.GetSubstructMatches(monosaccharide_pattern)) < 2:
        return False, "Insufficient monosaccharide units"

    # Check for glycosidic linkages connecting monosaccharides
    if len(mol.GetSubstructMatches(glycosidic_bond_pattern)) < 1:
        return False, "No glycosidic linkages found"
    
    return True, "Contains sufficient monosaccharide units and glycosidic linkages indicative of an oligosaccharide"

# This function aims to recognize the defining patterns of oligosaccharides; it is a basic implementation and might need further refinement for edge cases.
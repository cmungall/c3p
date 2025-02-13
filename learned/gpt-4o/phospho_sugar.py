"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar ring (5 or 6-membered ring with multiple hydroxyl groups)
    sugar_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1 | C1C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring structure found"
    
    # Look for phosphate ester linkage: -C-O-P(=O)(O)O
    phosphate_pattern = Chem.MolFromSmarts("[CX4,CX3]-[OX2]-P(=O)([OX1])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester linkage found"

    return True, "Contains sugar ring structure with a phosphate ester linkage"
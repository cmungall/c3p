"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:27600 phospho sugar
Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.

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
    
    # Check for monosaccharide backbone
    monosaccharide_pattern = Chem.MolFromSmarts("[OX2]([CX4]([CX4]([OX2])[CX3](=[OX1])[OX2])[CX3](=[OX1])[CX4]([CX4]([OX2])[OX2])[OX2])[H]")
    if not mol.HasSubstructMatch(monosaccharide_pattern):
        return False, "No monosaccharide backbone found"
    
    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Check for ester bond between phosphate and alcohol group
    ester_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])[OX2]")
    if not AllChem.FindMolChemismatch(mol, ester_pattern):
        return False, "No ester bond between phosphate and alcohol group"
    
    return True, "Contains a monosaccharide backbone with a phosphate group esterified to an alcoholic hydroxy group"
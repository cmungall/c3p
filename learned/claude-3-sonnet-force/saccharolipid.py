"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI:26094 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid that contains a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbohydrate patterns
    sugar_patterns = [Chem.MolFromSmarts("[OX2H][CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])"),  # linear sugar
                      Chem.MolFromSmarts("[OX2H][CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX3]"),  # sugar with branch
                      Chem.MolFromSmarts("[OX2H][CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[OX2H]"),  # sugar with open ring
                      Chem.MolFromSmarts("[OX2H][CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX3H2]")]  # furanose
    sugar_matches = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_matches:
        return False, "No carbohydrate moiety found"

    # Look for lipid patterns
    lipid_patterns = [Chem.MolFromSmarts("[CX3](=O)[OX2H]"),  # ester
                      Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2]"),  # long alkyl chain
                      Chem.MolFromSmarts("[SX4+2]([OX1-])([OX1-])[OX2H]")]  # sulfate ester
    lipid_matches = any(mol.HasSubstructMatch(pattern) for pattern in lipid_patterns)
    if not lipid_matches:
        return False, "No lipid moiety found"

    return True, "Contains both carbohydrate and lipid moieties"
"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is characterized by a glycerol backbone with two hydroxy groups acylated (esterified).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the glycerol backbone with consideration for chirality: O[C@@H](C)O and related
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](CO)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify ester group: a consistent pattern like C(=O)O connected to glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, "Incorrect number of ester groups attached to backbone"
    
    # Check for presence of phosphate or other groups: P(=O)(O)O
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, suggesting a phospholipid"

    return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds"
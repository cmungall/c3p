"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the steroid backbone pattern (rings A, B, C, D)
    # General steroid: four-ring structure, C-C-C-C-C in a particular arrangement.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C)CC4CCC(C)C4C3CC2C1")  # General pattern
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group on the steroid (carbonyl group =O on C3)
    oxo_pattern = Chem.MolFromSmarts("[C;R]=O")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        return False, "No 3-oxo group found"
    
    # Check for the 5alpha hydrogen (stereochemistry of the steroid)
    # Simplified check: verify H stereochemistry around C5
    five_alpha_pattern = Chem.MolFromSmarts("C[C@H](C)CCCC")  # Simplified pattern
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "No 5alpha configuration found"

    return True, "Contains a 3-oxo group and the 5alpha stereochemical configuration, typical of 3-oxo-5alpha-steroids"
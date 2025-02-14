"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a glycerol in which the acyl substituent is located
    specifically at the hydroxyl group of the second carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with acyl group at the second position
    glycerol_2_monoglyceride_pattern = Chem.MolFromSmarts("[CX4](CO)(CO)OC(=O)C")
    
    if mol.HasSubstructMatch(glycerol_2_monoglyceride_pattern):
        # Count ester groups, should be exactly one
        ester_pattern = Chem.MolFromSmarts("OC(=O)C")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) != 1:
            return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
        
        return True, "Contains a glycerol backbone with acyl substituent at position 2"
    
    return False, "Does not match the 2-monoglyceride pattern"
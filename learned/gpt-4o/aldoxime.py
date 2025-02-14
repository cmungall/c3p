"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is an oxime of an aldehyde RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldoxime pattern C=N-OH where C is double-bonded to N and single-bonded to H
    aldoxime_pattern = Chem.MolFromSmarts("[CX3;$([#6])](=[NX2][OX1H])")

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains C=N-OH group with C bonded to a hydrogen"

    return False, "Does not contain the C=N-OH group with C bonded to a hydrogen"

# Sample test
smiles = "C(CSC)=NO"  # example of an aldoxime
result, reason = is_aldoxime(smiles)
print(result, reason)
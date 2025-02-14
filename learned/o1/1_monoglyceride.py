"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI:17854 1-monoglyceride
"""

from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 1-monoglyceride
    # Pattern for glycerol backbone with ester at position 1
    # [C](=O)OC[CH](O)C(O)O
    pattern_smarts = '[C](=O)O[CH2][CH](O)[CH2]O'
    pattern = Chem.MolFromSmarts(pattern_smarts)

    if mol.HasSubstructMatch(pattern):
        # Check that there is only one ester group
        ester_smarts = '[CX3](=O)[OX2H0][#6]'
        ester_pattern = Chem.MolFromSmarts(ester_smarts)
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) == 1:
            return True, "Contains glycerol backbone with acyl group at position 1"
        else:
            return False, f"Found {len(ester_matches)} ester groups, requires exactly 1"
    else:
        return False, "Does not match 1-monoglyceride pattern"
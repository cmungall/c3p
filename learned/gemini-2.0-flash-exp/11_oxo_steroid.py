"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid core, focusing on ring structure and relative position of C11
    # The goal is to capture the ABCD ring system, and a carbonyl group on C11
    # The pattern is built stepwise to make it readable
    # Core of the steroid is a 6-6-6-5 fused ring system
    # [C](-[C]-) is to capture the C10-C19 bond, where the 19 methyl can be present or not
    # [C](-[C]-) is to capture the C13-C18 bond, where the 18 methyl can be present or not
    # [C](-[C]-) is to capture the C9-C11 bond, C9 should have two carbons attached to it, while C11 should have one carbonyl, and 2 other carbons attached to it
    # [C] is to capture C12
    # [C](-[C]-) is to capture C8-C14 bond
    # This is a simplification, but should capture the essential connectivity and the position 11 carbonyl
    steroid_core_smarts = "[C]12[C](-[C]-)[C]3[C](-[C]-)[C](=[O])[C](-[C]-)[C](-[C]-)[C]4[C]3[C]2[C](-[C]-)[C]14"

    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    if not mol.HasSubstructMatch(steroid_core_pattern):
          return False, "Not a steroid with carbonyl at position 11"

    
    return True, "Steroid with a carbonyl group at position 11"
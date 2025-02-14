"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid has a 3-ketone group and a C=C double bond between position 4 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for steroid core, 3-oxo group, and delta(4) double bond
    # The core is a fused ring system and the numbered carbons correspond to a steroid numbering scheme:
    #      1   2   3   4
    #    /---\---/-\---/
    #   10    9   8   5
    #   \---\---/-\---/
    #      6   7
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]([C]1)[C]([C])([C])[C]3[C]2[C][C][C]4[C]3([C])[C][C][C]4")
    oxo_group_pattern = Chem.MolFromSmarts("[C]1[C][C](=[O])[C]=[C]1")
    # the conjugation of the C=C double bond with the carbonyl is implicitly included in the oxo_group_pattern
    
    # Check for substructure matches
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo-Delta(4) group found"

    return True, "Contains steroid core with 3-oxo and Delta(4) double bond."
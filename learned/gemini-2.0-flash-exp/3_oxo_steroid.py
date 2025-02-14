"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid with a carbonyl group (=O) at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for steroid core (simplified, focuses on ring connectivity).
    # The key here is to identify a pattern that contains 3 fused six membered rings and one five membered ring
    # and at least one carbon with at least 2 neighbors
    steroid_core_pattern = Chem.MolFromSmarts("[CH2X4,CHX3,CHX2]1[CH2X4,CHX3,CHX2]2[CH2X4,CHX3,CHX2]3[CH2X4,CHX3,CHX2][CH2X4,CHX3,CHX2]4[CH2X4,CHX3,CHX2]1[CH2X4,CHX3,CHX2]5[CH2X4,CHX3,CHX2]6[CH2X4,CHX3,CHX2]2[CH2X4,CHX3,CHX2]7[CH2X4,CHX3,CHX2]3[CH2X4,CHX3,CHX2]4[CH2X4,CHX3,CHX2]756")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have a steroid core"

    # SMARTS for carbonyl at position 3. It will be attached to a carbon of the 6 membered ring
    # connected to the fused ring.
    # This focuses on a carbon being part of a 6 membered ring and having a double bonded oxygen.
    oxo_at_3_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4]([CH2X4,CHX3,CHX2])[CH2X4,CHX3,CHX2][CH2X4,CHX3,CHX2][CH2X4,CHX3,CHX2][CH2X4,CHX3,CHX2][CH2X4,CHX3,CHX2]")

    if not mol.HasSubstructMatch(oxo_at_3_pattern):
         return False, "Molecule does not have a carbonyl at position 3"
    
    return True, "Molecule is a 3-oxo steroid"
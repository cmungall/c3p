"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1)-steroid based on its SMILES string.
    A 3-oxo-Delta(1)-steroid is a steroid with a ketone at position 3 and a double bond
    between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1)-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the steroid core, 3-oxo group and delta(1) double bond.
    # The numbers in square brackets represent the ring ids to match exact position on steroid core.
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]3[C]([C]12)[C][C][C]4[C]3([C][C]4)")
    oxo_3_pattern = Chem.MolFromSmarts("[C]1(=[O])[C][C][C]2[C]3[C]([C]12)[C][C][C]4[C]3([C][C]4)")
    delta_1_pattern = Chem.MolFromSmarts("[C]1=C[C][C]2[C]3[C]([C]12)[C][C][C]4[C]3([C][C]4)")
    delta_1_oxo_3_pattern = Chem.MolFromSmarts("[C]1(=[O])=[C][C][C]2[C]3[C]([C]12)[C][C][C]4[C]3([C][C]4)")
    
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found."
    if not mol.HasSubstructMatch(oxo_3_pattern):
         return False, "No 3-oxo group found in the steroid core."
    if not mol.HasSubstructMatch(delta_1_pattern):
        return False, "No double bond between positions 1 and 2 in steroid core."
    if not mol.HasSubstructMatch(delta_1_oxo_3_pattern):
          return False, "Molecule is not a 3-oxo-Delta(1)-steroid, double bond and oxo group not in correct position"    

    return True, "Molecule is a 3-oxo-Delta(1)-steroid."
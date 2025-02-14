"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid core (simplified)
    # The exact positions are handled in the other SMARTS queries
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C]([C])[C]3[C]2[C][C][C]4[C]3[C][C]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core not found"

    # Define SMARTS pattern for the 3-beta-hydroxyl group
    # We use [C@H] to specify the stereo configuration of the 3 position hydrogen
    # which needs to be on the opposite side of the beta hydroxyl group.
    beta_hydroxy_pattern = Chem.MolFromSmarts("[C@H]1[C]([OH])[C]([C])[C]2[C]1[C]([C])[C]3[C]2[C][C][C]4[C]3[C][C]4")

    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
          return False, "No 3-beta hydroxyl group found"

    # Define SMARTS pattern for the double bond between C5 and C6
    # Use [C]=C
    delta5_bond_pattern = Chem.MolFromSmarts("[C]1[C]([C]=[C])[C]([C])[C]2[C]1[C]([C])[C]3[C]2[C][C][C]4[C]3[C][C]4")
    if not mol.HasSubstructMatch(delta5_bond_pattern):
         return False, "No double bond between C5 and C6"

    return True, "3beta-hydroxy-Delta(5)-steroid identified"
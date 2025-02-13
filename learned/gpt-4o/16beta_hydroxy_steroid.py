"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is a hydroxy steroid with an -OH group at position 16 in beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid backbone check (e.g., four rings fused)
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C4(C)CCC(O)CC4CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for beta-OH at position 16
    beta_oh_16_pattern = Chem.MolFromSmarts("[C@H](O)[CH](C)[CH2][CH2][CH2][CH2]C")
    if not mol.HasSubstructMatch(beta_oh_16_pattern):
        return False, "No beta-OH at position 16 found"

    return True, "Contains beta-hydroxy group on steroid nucleus at position 16"
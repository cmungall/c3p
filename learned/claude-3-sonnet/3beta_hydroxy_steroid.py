"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36621 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern
    # Steroid backbone is a fused ring system with 3 6-membered rings and 1 5-membered ring
    steroid_pattern = Chem.MolFromSmarts('[C@@]1(CC[C@@]2([H])[C@]3([H])CC[C@@]4([H])C[C@@]([H])(CC[C@@]4([H])C3)C2)CC[C@]1([H])O')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl group at position 3 in beta configuration
    hydroxy_pattern = Chem.MolFromSmarts('[C@@]1(CC[C@@]2([H])[C@]3([H])CC[C@@]4([H])C[C@@]([H])(CC[C@@]4([H])C3)[C@H]2O)CC[C@]1([H])O')
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxyl group found"

    return True, "Contains a steroid backbone with a 3beta-hydroxyl group"
"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is a 3beta-hydroxy-steroid that contains a double bond
    between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 3beta-hydroxy group SMARTS pattern
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]1([OH])[C@@H]2CC[C@]3(CC[C@@H]4CCC=C4C3)[C@@H]2C[C@]1(C)C")

    # Steroid backbone (3 six-membered rings and 1 five-membered ring)
    steroid_pattern = Chem.MolFromSmarts(
        "[C&R2]12[C@@H](C[C@H]3CC[C@@H]4CC[C@@H](O)CC4=C3C1)C=C2"
    )

    # Delta(5) double bond pattern
    delta5_pattern = Chem.MolFromSmarts("C=C")

    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for Delta(5) double bond
    if not mol.HasSubstructMatch(delta5_pattern):
        return False, "No Delta(5) double bond found"

    return True, "Molecule is a 3beta-hydroxy-Delta(5)-steroid"
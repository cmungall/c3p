"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is characterized by a hydroxysteroid backbone with attached sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for steroid and glycoside features
    steroid_pattern = Chem.MolFromSmarts("C1CC2=C3C(C=C4C3CCC4)CC2C1")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    glycoside_pattern = Chem.MolFromSmarts("[O;D2]-;!@[C;R]")  # Sugars typically attached via an oxygen atom

    # Check for steroid backbone with hydroxyl groups
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups on steroid backbone"

    # Check for attached glycoside(s)
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycoside linkages found"

    # If all patterns are matched, classify as steroid saponin
    return True, "Contains hydroxysteroid backbone with glycoside linkages"
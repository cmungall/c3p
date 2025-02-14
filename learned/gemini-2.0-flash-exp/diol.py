"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a molecule with exactly two hydroxyl (-OH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Find hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]") # oxygen with 2 bonds and 1 hydrogen

    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    num_hydroxyls = len(hydroxyl_matches)

    # Check for exactly two hydroxyl groups
    if num_hydroxyls == 2:
        return True, "Molecule contains exactly two hydroxyl groups"
    elif num_hydroxyls < 2:
        return False, f"Molecule has fewer than two hydroxyl groups, it has {num_hydroxyls}."
    else:
        return False, f"Molecule has more than two hydroxyl groups, it has {num_hydroxyls}."
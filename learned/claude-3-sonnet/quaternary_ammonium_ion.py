"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35485 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules containing phosphate or phosphocholine groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("OP(O)([O-])O")):
        return False, "Contains phosphate or phosphocholine group"

    # Check for quaternary nitrogen with four organyl substituents
    quat_n_pattern = Chem.MolFromSmarts("[N+]([C,c])([C,c])([C,c])([C,c])")
    if mol.HasSubstructMatch(quat_n_pattern):
        return True, "Contains quaternary nitrogen atom with four organyl substituents"

    return False, "No quaternary nitrogen atom with four organyl substituents found"
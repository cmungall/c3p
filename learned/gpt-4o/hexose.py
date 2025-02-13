"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is defined as a six-carbon monosaccharide with an aldehyde at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count all carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Ensure exactly six carbon atoms are present
    if c_count != 6:
        return False, f"Contains {c_count} carbon atoms, but a hexose requires exactly 6"

    # Define SMARTS patterns for functional groups in linear and cyclic form
    aldehyde_pattern = Chem.MolFromSmarts("O=CC(C)(C)C")
    ketone_pattern = Chem.MolFromSmarts("C(C)(C)C(=O)C")

    # Check for linear forms: aldehyde or ketone groups
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Structure matches an aldohexose (aldehyde group at position 1)"
    
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Structure matches a ketohexose (ketone group at position 2)"

    # Enhanced patterns for furanose and pyranose
    furanose_pattern = Chem.MolFromSmarts("[C@H]1O[C@@H](C(C1)O)O")
    pyranose_pattern = Chem.MolFromSmarts("[C@H]1O[C@@H](C(O)[C@H](O)C1)O")

    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        return True, "Structure matches a cyclic furanose or pyranose form"

    return False, "Does not contain hexose-defining functional groups (aldehyde or ketone)"
"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a potential aldehyde group or existing in
    a cyclic (furanose or pyranose) form characteristic of aldopentoses.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"

    # Detect potential linear aldopentose with an aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](=O)[C;$(C[OH])]")  # more specific terminal ald group pattern
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains 5 carbon atoms and a terminal aldehyde group"

    # Detect potential pyranose form (six-membered ring with 5 carbons and 1 oxygen)
    pyranose_pattern = Chem.MolFromSmarts("C1C[C@H](O)[C@H](O)[C@H](O)O1")
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains 5 carbon atoms forming a typical pyranose ring"

    # Detect potential furanose form (five-membered ring with 4 carbons and 1 oxygen)
    furanose_pattern = Chem.MolFromSmarts("C1C[C@@H](O)[C@@H](O)O1")
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains 5 carbon atoms forming a typical furanose ring"

    return False, "Does not match an aldopentose structure"
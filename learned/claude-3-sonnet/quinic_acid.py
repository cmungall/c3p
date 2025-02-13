"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:16857 quinic acid
A cyclitol carboxylic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid, with a cyclohexane core,
    multiple hydroxyl groups, and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclohexane core
    cyclohexane_pattern = Chem.MolFromSmarts("[C&R1]1[C&R1][C&R1][C&R1][C&R1][C&R1]1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane core found"

    # Look for at least three hydroxyl groups (cyclitol)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Fewer than 3 hydroxyl groups found"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Oxygen count too low for quinic acid"

    return True, "Contains cyclohexane core with multiple hydroxyl groups and a carboxylic acid"
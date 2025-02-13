"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: CHEBI:78354 anthoxanthin

Anthoxanthins are a type of flavonoid pigments in plants. They are water-soluble pigments that range in color from white or colorless to a creamy yellow, often found in petals of flowers.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavone backbone
    flavone_pattern = Chem.MolFromSmarts("c1cc(=O)c2c(cccc2)o1")
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No flavone backbone found"

    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Fewer than 2 hydroxyl groups found"

    # Look for glycosidic bonds (oxygen attached to two carbons)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2;!R]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    # Check for other substituents like methoxy, sulfate, etc.
    substituents = ["[OX2C]", "[SX4](=O)(=O)([OX1])", "[CX3](=O)[OX2H]"]
    sub_matches = []
    for sub_pattern in substituents:
        sub_pattern = Chem.MolFromSmarts(sub_pattern)
        sub_matches.extend(mol.GetSubstructMatches(sub_pattern))

    # Anthoxanthins are typically large, colored, and water-soluble
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 300:
        return False, "Molecular weight too low for anthoxanthin"

    if len(glycosidic_matches) > 0 or len(sub_matches) > 0:
        return True, "Contains flavone backbone with hydroxyl and glycosidic/other substituents"
    else:
        return False, "No glycosidic or other common substituents found"
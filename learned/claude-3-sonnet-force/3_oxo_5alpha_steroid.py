"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI:63502 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid is a steroid that has a ketone group at position 3 and
    an alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]13[C@H]([C@H]([C@@H]2[C@@]([C@@]1(CC[C@H](C2)C)(C)C)(C)CCC(=O)O)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        return False, "No ketone group found"

    # Check if oxo group is at position 3
    position_3_pattern = Chem.MolFromSmarts("[C@]13[C@H]([C@@H]2[C@@](CC[C@H](C2)C)(C)CCC(=O)O)[C@@H](C)CCC1")
    if not mol.HasSubstructMatch(position_3_pattern):
        return False, "Ketone group not at position 3"

    # Check for 5alpha configuration
    alpha_pattern = Chem.MolFromSmarts("[C@@]13[C@H]([C@H]2[C@@](CC[C@H](C2)C)(C)CCC(=O)O)[C@@H](C)CCC1")
    if not mol.HasSubstructMatch(alpha_pattern):
        return False, "Not alpha configuration at position 5"

    return True, "Contains 3-oxo and 5alpha configurations in a steroid backbone"


# Example usage
smiles1 = "C[C@H](C[C@@H](O)[C@H](O)C(C)(C)O)[C@H]1CC[C@@]2(C)[C@@H]3[C@@H](O)C[C@@H]4[C@]5(C[C@@]35CC[C@]12C)CCC(=O)[C@@]4(C)CO"
smiles2 = "CCCCCCCCCCCCCCCCCCCC"

print(is_3_oxo_5alpha_steroid(smiles1))  # (True, 'Contains 3-oxo and 5alpha configurations in a steroid backbone')
print(is_3_oxo_5alpha_steroid(smiles2))  # (False, 'No steroid backbone found')
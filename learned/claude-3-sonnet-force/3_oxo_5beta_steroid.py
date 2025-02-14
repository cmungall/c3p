"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:17971 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ketone at position 3
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone group found at position 3"

    # Look for beta configuration at position 5
    beta_pattern = Chem.MolFromSmarts("[C@@H]1[C@H](CC[C@@]2([C@]1([C@H]([C@@]3([C@H](CC2)C)C)(C)C)C)C")
    if not mol.HasSubstructMatch(beta_pattern):
        return False, "No beta configuration found at position 5"

    # Check if molecule is a steroid
    steroid_pattern = Chem.MolFromSmarts("[C@]12CC[C@@]3([C@@]1(CCC[C@@H]2O)C)[C@H](CC[C@@H]4[C@]3(CCC(=O)C[C@@H]4)C)C"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Not a steroid backbone"

    return True, "Molecule contains a ketone at position 3 and beta configuration at position 5 on a steroid backbone"
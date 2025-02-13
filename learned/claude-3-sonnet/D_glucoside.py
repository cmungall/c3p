"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:18237 D-glucoside
A D-glucoside is any glucoside in which the glycoside group is derived from D-glucose.
"""

from rdkit import Chem
from rdkit.Chem import rdqueries, AllChem

def is_D_glucoside(smiles):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for D-glucose pattern
    d_glucose_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@@H]([C@H](C(O)=O)O)O)O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(d_glucose_pattern):
        return False, "No D-glucose substructure found"

    # Check for glycosidic bond (-O-C-O-) with the glucose
    glycosidic_bond_pattern = Chem.MolFromSmarts("[O;r]")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)

    glucose_atoms = list(map(lambda x: x[0], mol.GetSubstructMatches(d_glucose_pattern)))
    glycosidic_glucose_bonds = [bond for bond in glycosidic_bonds if bond in glucose_atoms]

    if not glycosidic_glucose_bonds:
        return False, "No glycosidic bonds found with D-glucose"

    return True, "Contains D-glucose substructure with a glycosidic bond"
"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:26722 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is defined as any 3-oxo steroid that contains a double bond between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@@]3([C@@]4([C@@]1(CC[C@@]2([H])C)[H])CC[C@H]5[C@@]4(C)CC[C@](C)(C5)C)[C@@]34C")
    if not mol.HasSubstructMatch(steroid_pattern, useChirality=True):
        return False, "No steroid backbone found"

    # Check for 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("[C](=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) != 1:
        return False, "Incorrect number of oxo groups"

    # Check if the oxo group is at the 3-position
    oxo_atom = mol.GetAtomWithIdx(list(oxo_matches[0])[0])
    ring_atoms = mol.GetRingAtoms()
    if oxo_atom not in ring_atoms:
        return False, "Oxo group not at the 3-position"

    # Check for Delta(1) double bond
    double_bond_pattern = Chem.MolFromSmarts("[C]=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 1:
        return False, "Incorrect number of double bonds"

    # Check if double bond is at Delta(1) position
    double_bond_atoms = [mol.GetAtomWithIdx(idx) for idx in list(double_bond_matches[0])]
    ring_atoms = mol.GetRingAtoms()
    if not all(atom in ring_atoms for atom in double_bond_atoms):
        return False, "Double bond not at Delta(1) position"

    return True, "Contains a 3-oxo group and a double bond at the Delta(1) position of the steroid backbone"
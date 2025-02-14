"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:37504 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondDir

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is an oxo steroid where an oxo substituent is located at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@@]1(CC[C@]2([H])[C@]3([H])CC[C@@H]4[C@@]([H])(CC[C@@]5([H])[C@@H](CC[C@@]45[H])C3)C2)[C@@H](C)CCC1"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for oxo group at position 3
    oxo_pattern = Chem.MolFromSmarts("[C@H](=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No oxo group found"

    # Get atom index of oxo group
    oxo_atom_idx = next((atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8), None)

    # Check if oxo group is attached to the right carbon (position 3)
    if oxo_atom_idx is None:
        return False, "No oxo group found"

    oxo_neighbor_idx = mol.GetAtomWithIdx(oxo_atom_idx).GetNeighbors()[0].GetIdx()
    oxo_neighbor_bond = mol.GetBondBetweenAtoms(oxo_atom_idx, oxo_neighbor_idx)
    if oxo_neighbor_bond.GetBondDir() != BondDir.ENDDOWNRIGHT:
        return False, "Oxo group is not at position 3"

    return True, "Contains an oxo group at position 3 of the steroid backbone"
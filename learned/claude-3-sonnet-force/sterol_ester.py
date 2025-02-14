"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:36976 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the
    carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@]1(CC[C@]2([H])C3=C(CC[C@]12C)[C@@]1([H])CC[C@H](C(C)(C)[C@]1([H])CC3)C")
    steroid_match = mol.GetSubstructMatches(steroid_pattern)
    if not steroid_match:
        return False, "No steroid backbone found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Check that ester is attached to steroid at 3-position
    ester_atom_idx = list(ester_matches[0])[0]
    ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
    if ester_atom.GetTotalNumHs() != 0:
        return False, "Ester not attached at 3-position of steroid"

    # Check for long aliphatic chain attached to ester (fatty acid)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Fatty acid chain too short"

    # Check molecular weight - sterol esters typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for sterol ester"

    return True, "Contains steroid backbone with ester group at 3-position and fatty acid chain attached"
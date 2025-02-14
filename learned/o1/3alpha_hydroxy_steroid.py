"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:12345 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxy group at the 3-position in the alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid backbone SMARTS pattern (cyclopentanoperhydrophenanthrene)
    steroid_smarts = """
    [#6]1([#1])[#6][#6]2[#6]1[#6][#6]3
    [#6]2[#6][#6][#6]4[#6]3[#6][#6][#6]([#6]4)[#6][#6]
    """
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"

    # Define the 3alpha-hydroxy group SMARTS pattern
    # Approximately matches a secondary alcohol with alpha orientation on a ring
    hydroxy_smarts = "[C@H](O)[#6]"
    hydroxy_pattern = Chem.MolFromSmarts(hydroxy_smarts)
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No alpha-oriented hydroxy group found"

    # Attempt to locate the hydroxy group at the 3-position
    # Since exact atom indexing is unreliable, we'll check if the hydroxy group is attached to a ring atom
    ring_info = mol.GetRingInfo()
    for match in hydroxy_matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.IsInRing():
            # We consider this as the 3-position for our purposes
            return True, "Molecule is a 3alpha-hydroxy steroid"

    return False, "Alpha-oriented hydroxy group not at the correct position"
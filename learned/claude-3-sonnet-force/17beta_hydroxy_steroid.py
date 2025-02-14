"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: CHEBI:35577 17beta-hydroxy steroid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Bond

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is a steroid with a hydroxy group at position 17 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]1(C)[C@H]2[C@@]3([C@H](C[C@]4([C@@]3(C[C@H](C2)C4(C)C)C)C)C)CC[C@@]1(O)C(=O)CO"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Find the 17-position oxygen (attached to a cyclohexane ring)
    o_atom = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            neighbors = [mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType() for nbr in atom.GetNeighbors()]
            if Bond.SINGLE in neighbors and Bond.SINGLE in neighbors:  # Connected to two carbon atoms
                ring_info = mol.GetRingInfo().IsBondInRingOfSize(mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[0].GetIdx()).GetIdx(), 6)
                if ring_info:  # One of the carbons is part of a cyclohexane ring
                    o_atom = atom
                    break

    if o_atom is None:
        return False, "No 17-hydroxy group found"

    # Check configuration of 17-hydroxy group
    conf = mol.GetConformer(-1)
    o_neighbors = [nbr for nbr in o_atom.GetNeighbors()]
    c1_idx = o_neighbors[0].GetIdx()
    c2_idx = o_neighbors[1].GetIdx()
    c1_pos = conf.GetAtomPosition(c1_idx)
    c2_pos = conf.GetAtomPosition(c2_idx)
    o_pos = conf.GetAtomPosition(o_atom.GetIdx())
    v1 = c1_pos - o_pos
    v2 = c2_pos - o_pos
    cross = v1.CrossProduct(v2)

    # Positive cross product indicates beta configuration
    if cross.z > 0:
        return True, "Molecule has a 17beta-hydroxy group on the steroid backbone"
    else:
        return False, "17-hydroxy group is in the alpha configuration"
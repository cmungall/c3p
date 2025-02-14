"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:34545 azole
An azole is any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
Azoles can also contain one or more other non-carbon atoms, such as nitrogen, sulfur or oxygen.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a 5-membered ring containing nitrogen
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    azole_rings = []
    for ring in rings:
        if len(ring) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            azole_rings.append(ring)

    if not azole_rings:
        return False, "No 5-membered ring containing nitrogen found"

    # Check if the ring is aromatic
    for ring in azole_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_bond_types = [mol.GetBondBetweenAtoms(ring_atoms[i].GetIdx(), ring_atoms[i+1].GetIdx()).GetBondType()
                           for i in range(len(ring_atoms)-1)]
        ring_bond_types.append(mol.GetBondBetweenAtoms(ring_atoms[-1].GetIdx(), ring_atoms[0].GetIdx()).GetBondType())
        if all(bond_type == Chem.BondType.AROMATIC for bond_type in ring_bond_types):
            break
    else:
        return False, "5-membered nitrogen ring is not aromatic"

    # Check for allowed heteroatoms (N, O, S)
    allowed_heteroatoms = [7, 8, 16]
    ring_heteroatoms = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() not in [1, 6]]
    if any(heteroatom not in allowed_heteroatoms for heteroatom in ring_heteroatoms):
        return False, "Ring contains disallowed heteroatoms"

    return True, "Molecule contains a monocyclic aromatic 5-membered ring with nitrogen and possibly other allowed heteroatoms (O, S)"
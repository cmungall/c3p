"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:24636 azole
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

    # Look for 5-membered rings containing nitrogen
    azole_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if len(ring) == 5 and any(atom.GetAtomicNum() == 7 for atom in atoms):
            azole_rings.append(ring)

    if not azole_rings:
        return False, "No 5-membered rings containing nitrogen found"

    # Check if any of the 5-membered rings are aromatic
    is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for ring in azole_rings for idx in ring)
    if not is_aromatic:
        return False, "No aromatic 5-membered rings containing nitrogen found"

    # Check for other heteroatoms (N, O, S)
    heteroatoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    if not any(x in heteroatoms for x in [7, 8, 16]):
        return False, "No additional heteroatoms (N, O, S) found"

    return True, "Molecule contains an aromatic 5-membered ring with nitrogen and other heteroatoms"
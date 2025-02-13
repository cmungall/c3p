"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: CHEBI:51584 tetrahydrofuranone
Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for tetrahydrofuran ring
    thf_pattern = Chem.MolFromSmarts("C1CCOC1")
    thf_match = mol.GetSubstructMatches(thf_pattern)
    if not thf_match:
        return False, "No tetrahydrofuran ring found"

    # Check if there is an oxo substituent on the ring
    oxo_pattern = Chem.MolFromSmarts("C1CCOC(=O)1")
    oxo_match = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_match:
        return False, "No oxo substituent on the tetrahydrofuran ring"

    # Additional checks
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return False, "More than one ring present"

    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != 5:
        return False, "Ring size is not 5"

    # Count carbon and oxygen atoms in the ring
    c_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
    if c_count != 4 or o_count != 1:
        return False, "Ring does not contain 4 carbons and 1 oxygen"

    return True, "Molecule is a tetrahydrofuranone"
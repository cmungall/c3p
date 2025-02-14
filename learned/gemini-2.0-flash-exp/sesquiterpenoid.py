"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a 15-carbon skeleton derived from three isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for 15 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, f"Incorrect number of carbons: {c_count} (should be 15)"

    # 2. Check for isoprene units (using SMARTS pattern).  Allow carbons to be CH2 or CH or C.
    #   a basic branched carbon pattern can be used to detect isoprenoid structures.
    isoprene_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])([CX4])")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:  # Typically expect at least 2 branched centers
            return False, f"Not branched as expected for isoprene units, only {len(isoprene_matches)} found"

    # 3. Check for cyclic structures (at least one ring).
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No ring structure found"


    return True, "Contains 15 carbons, appears to be built from isoprene units, and contains a ring."
"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid typically has a core terpenoid structure with 20 carbons,
    potentially rearranged or functionalized.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Insufficient carbon atoms for diterpenoid (less than 20)"

    # Family specific patterns are typically challenging, hence our focus will be on
    # ring structures and diverse functionalities typical for diterpenoids.
    # Check for the presence of complex ring structures
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or ring_info.NumRings() == 0:
        return False, "No ring structures found"

    # Check for specific terpenoid-like ring patterns
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size >= 6 for size in ring_sizes):
        return False, "Lacks larger ring structures typical of diterpenoids"

    # Additional checks can include specific functional groups or structural motifs
    # Present in known diterpenoids, though complex variations exist.

    # We consider here the generic properties without a highly specific pattern due to diversity
    return True, "Matches general diterpenoid structure based on ring presence and carbon count"
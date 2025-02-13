"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic aromatic hydrocarbon (PAH) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate number of aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = [ring for ring in ring_info.BondRings() if all(mol.GetBondBetweenAtoms(a1, a2).GetIsAromatic() for a1, a2 in zip(ring, ring[1:] + ring[:1]))]

    # Check if there is more than one aromatic ring
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings, not a polycyclic arene"

    # Check for presence of fused aromatic rings (shared edges)
    has_fused_rings = any(len(set(aromatic_rings[i]).intersection(aromatic_rings[j])) > 1 for i in range(len(aromatic_rings)) for j in range(i+1, len(aromatic_rings)))
    
    if not has_fused_rings:
        return False, "No fused aromatic rings, not a polycyclic arene"
    
    return True, "Structure is consistent with a polycyclic aromatic hydrocarbon"
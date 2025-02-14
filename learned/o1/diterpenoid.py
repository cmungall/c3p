"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:35196 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is derived from a diterpene (composed of four isoprene units),
    but may have rearranged or modified skeletons, often missing some methyl groups.

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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 25:
        return False, f"Carbon count ({c_count}) not in typical range for diterpenoids (15-25 carbons)"
    
    # Attempt to identify isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, less than expected for diterpenoids"

    # Check for terpenoid-like methyl branching
    methyl_branch_pattern = Chem.MolFromSmarts("C(C)(C)C")
    methyl_branch_matches = mol.GetSubstructMatches(methyl_branch_pattern)
    if len(methyl_branch_matches) == 0:
        return False, "No methyl branching found typical of terpenoids"

    # Consider rings (many diterpenoids are cyclic)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No rings found; diterpenoids often contain ring structures"

    # Check for oxygen-containing functional groups (common in diterpenoids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; diterpenoids are often oxidized"

    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.2f}) not in typical range for diterpenoids (250-500 Da)"

    return True, "Molecule has features consistent with diterpenoids (carbon count, isoprene units, methyl branching, rings, oxygen atoms)"
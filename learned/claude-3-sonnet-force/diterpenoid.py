"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:35924 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is any terpenoid derived from a diterpene, which may have undergone
    rearrangements or modifications of the C20 skeleton, such as the removal of methyl groups.

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

    # Check molecular weight range (typically 280-360 Da for diterpenoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 280 or mol_wt > 360:
        return False, f"Molecular weight {mol_wt:.2f} outside typical diterpenoid range"

    # Check for diterpenoid backbone (C20 skeleton)
    diterpenoid_backbone = Chem.MolFromSmarts("[C&r5,r6,r7,r8]")
    if not mol.HasSubstructMatch(diterpenoid_backbone):
        return False, "No diterpenoid backbone found"

    # Look for typical functional groups (hydroxy, carbonyl, ether, ester, etc.)
    functional_groups = Chem.MolFromSmarts("[OH,O,C(=O)O,COC]")
    if not mol.HasSubstructMatch(functional_groups):
        return False, "No typical diterpenoid functional groups found"

    # Check for ring systems (bicyclic or tricyclic)
    ring_info = mol.GetRingInfo()
    if not any(len(ring) >= 5 for ring in ring_info.BondRings()):
        return False, "No ring systems found"

    # Count carbon, hydrogen, and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 20 or c_count > 25:
        return False, f"Carbon count {c_count} outside typical diterpenoid range"
    if h_count < 24 or h_count > 40:
        return False, f"Hydrogen count {h_count} outside typical diterpenoid range"
    if o_count < 1 or o_count > 5:
        return False, f"Oxygen count {o_count} outside typical diterpenoid range"

    return True, "Molecule matches diterpenoid structural features"
"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:37242 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is any terpenoid derived from a diterpene (C20 skeleton),
    which may be rearranged or modified by the removal of methyl groups.

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

    # Check for C20 skeleton
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:
        return False, "Molecule has less than 20 atoms, cannot be a diterpenoid"

    # Check for presence of only C, H, and O atoms
    allowed_atoms = {6, 1, 8}  # C, H, O
    atoms = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if not atoms.issubset(allowed_atoms):
        return False, "Molecule contains atoms other than C, H, O"

    # Count rings
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings < 4:
        return False, "Diterpenoids typically have 4 or more rings"

    # Count rotatable bonds
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds < 2:
        return False, "Diterpenoids typically have 2 or more rotatable bonds"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 400:
        return False, "Molecular weight outside typical range for diterpenoids"

    return True, "Molecule meets structural requirements for a diterpenoid"
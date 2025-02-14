"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33587 polycyclic arene
A polycyclic aromatic hydrocarbon.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon.

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

    # Check for aromaticity
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    num_aromatic_heterocycles = rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
    if num_aromatic_rings == 0 and num_aromatic_heterocycles == 0:
        return False, "Molecule does not contain any aromatic rings"

    # Check for polycyclic structure
    fused_rings = mol.GetRingInfo().BondRings()
    if len(fused_rings) < 2:
        return False, "Molecule does not have fused rings"

    # Check for characteristic ring sizes and fused ring systems
    ring_sizes = [len(ring) for ring in fused_rings]
    if not any(size in [5, 6, 7] for size in ring_sizes):
        return False, "Molecule does not contain characteristic ring sizes for polycyclic arenes"

    # Check for only carbon and hydrogen atoms (optional)
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if any(atom_num not in [1, 6] for atom_num in atoms):
        return True, "Polycyclic aromatic hydrocarbon containing heteroatoms"

    return True, "Polycyclic aromatic hydrocarbon"
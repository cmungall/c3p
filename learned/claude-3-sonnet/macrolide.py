"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: CHEBI:35608 macrolide
A macrocyclic lactone with a ring of twelve or more members derived from a polyketide.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get largest ring size
    ring_info = mol.GetRingInfo()
    largest_ring_size = max([ring.AtomCount() for ring in ring_info.AtomRings()])

    # Check if molecule has a macrocyclic ring (>= 12 atoms)
    if largest_ring_size < 12:
        return False, "No macrocyclic ring found (largest ring size < 12)"

    # Check if molecule contains a lactone group
    lactone_pattern = Chem.MolFromSmarts("[C@H]1[C@@]([C@@H](C(=O)O1))(O)")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"

    # Check if molecule is derived from polyketide
    # Estimate based on presence of long carbon chains and carbonyl groups
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_carbonyl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(bond.GetBondTypeAsDouble() == 2 for bond in atom.GetBonds()) == 1)

    if n_carbons < 15 or n_oxygens < 4 or n_carbonyl < 2:
        return False, "Structure does not appear to be derived from polyketide"

    return True, "Contains macrocyclic lactone ring (>= 12 atoms) derived from polyketide"
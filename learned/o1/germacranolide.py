"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 15 carbons (sesquiterpene)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, f"Number of carbons is {c_count}, not 15 (not a sesquiterpene)"

    # Check for 10-membered ring (germacrane skeleton)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if 10 not in ring_sizes:
        return False, "No 10-membered ring found (germacrane skeleton missing)"

    # Check for lactone group (cyclic ester within a ring)
    lactone_smarts = "[C;R](=O)[O;R][C;R]"  # Ester group where atoms are in a ring
    lactone_mol = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_mol):
        return False, "No lactone ring (cyclic ester) found"

    return True, "Molecule is a germacranolide (15 carbons, 10-membered ring, and lactone group)"
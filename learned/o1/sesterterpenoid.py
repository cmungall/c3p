"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene (C25 backbone),
    possibly modified by rearrangement or removal of small groups (e.g., methyl groups).
    They typically contain multiple rings, methyl branching, and oxygen-containing functional groups.
    They are composed mainly of carbon and hydrogen, with oxygen as the main heteroatom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for heteroatoms other than oxygen
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 6 and atomic_num != 1 and atomic_num != 8:
            return False, "Contains heteroatoms other than carbon, hydrogen, or oxygen"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 22 or c_count > 28:
        return False, f"Carbon count is {c_count}, not typical for sesterterpenoids (expected ~25 carbons)"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; sesterterpenoids often contain oxygenated functional groups"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No ring structures found; sesterterpenoids often contain multiple rings"

    # Check for methyl groups (CH3)
    methyl_pattern = Chem.MolFromSmarts('[CX4H3]')
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 3:
        return False, f"Found {len(methyl_matches)} methyl groups; sesterterpenoids typically have multiple methyl branches"

    # Check for isopropyl groups (common in terpenoids)
    isopropyl_pattern = Chem.MolFromSmarts('[CH](C)C')
    isopropyl_matches = mol.GetSubstructMatches(isopropyl_pattern)
    if len(isopropyl_matches) < 1:
        return False, "No isopropyl groups found; terpenoids often contain isopropyl groups"

    # Check for terpenoid core using Murcko scaffold
    from rdkit.Chem.Scaffolds import MurckoScaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold)
    if scaffold.GetNumAtoms() < 10:
        return False, "Scaffold too small to be a terpenoid core"

    # Check for common terpenoid functional groups
    functional_groups = [
        Chem.MolFromSmarts('[CX3]=[OX1]'),    # Carbonyl group
        Chem.MolFromSmarts('[CX3](=O)[OX2H1]'),  # Carboxylic acid
        Chem.MolFromSmarts('[OX2H]'),         # Hydroxyl group
        Chem.MolFromSmarts('[CX3](=O)[OX2][CX3]'), # Ester
        Chem.MolFromSmarts('[OX2][CX4][OX2]'),    # Ether
    ]
    fg_found = False
    for fg in functional_groups:
        if mol.HasSubstructMatch(fg):
            fg_found = True
            break
    if not fg_found:
        return False, "No typical terpenoid functional groups found"

    return True, "Molecule meets criteria for a sesterterpenoid"
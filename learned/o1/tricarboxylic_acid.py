"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing exactly three carboxy groups,
    without peptide bonds or other interfering functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxy group pattern (both protonated and deprotonated forms)
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[O;H1,-1]')
    if carboxy_pattern is None:
        return False, "Failed to create carboxy group pattern"

    # Find all carboxy groups in the molecule
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    num_carboxy_groups = len(carboxy_matches)

    if num_carboxy_groups != 3:
        return False, f"Found {num_carboxy_groups} carboxy groups, need exactly 3"

    # Check for peptide bonds (amide linkages)
    peptide_bond_pattern = Chem.MolFromSmarts('[NX3][CX3](=O)[#6]')
    has_peptide_bond = mol.HasSubstructMatch(peptide_bond_pattern)
    if has_peptide_bond:
        return False, "Contains peptide bonds (amide linkages), not a simple tricarboxylic acid"

    # Check for ester groups (to exclude carboxylic acid derivatives)
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if ester_matches:
        return False, "Contains ester groups, carboxy groups are not free"

    # Check for size of the molecule (exclude very large molecules)
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 50:
        return False, f"Molecule is too large ({num_atoms} atoms), likely not a simple tricarboxylic acid"

    # Ensure carboxy groups are not part of ring systems (optional, depending on definition)
    for match in carboxy_matches:
        carbon_index = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        if carbon_atom.IsInRing():
            return False, "Carboxy group is part of a ring system"

    return True, "Contains exactly three free carboxy groups and no peptide bonds (tricarboxylic acid)"
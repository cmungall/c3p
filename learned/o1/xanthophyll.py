"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: xanthophyll
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygen-containing carotenoid; that is, a carotene (C40 polyene hydrocarbon) with oxygen-containing functional groups such as hydroxyls, ketones, epoxides, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbons ({c_count}) to be a carotenoid"

    # Check for extensive conjugation (long polyene chain)
    # Estimate the number of conjugated double bonds
    double_bond_count = 0
    conjugated_double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            # Check if connected to another double bond (conjugation)
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            for nbr in begin_atom.GetNeighbors():
                if nbr.GetIdx() != end_atom.GetIdx():
                    nbr_bond = mol.GetBondBetweenAtoms(begin_atom.GetIdx(), nbr.GetIdx())
                    if nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        conjugated_double_bond_count += 1
            for nbr in end_atom.GetNeighbors():
                if nbr.GetIdx() != begin_atom.GetIdx():
                    nbr_bond = mol.GetBondBetweenAtoms(end_atom.GetIdx(), nbr.GetIdx())
                    if nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        conjugated_double_bond_count += 1

    if conjugated_double_bond_count < 7:
        return False, f"Too few conjugated double bonds ({conjugated_double_bond_count}) for a carotenoid"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Check for typical oxygen-containing functional groups
    # Hydroxyl groups (-OH)
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]'))
    # Ketone groups (>C=O)
    has_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[#6]'))
    # Epoxide groups (C1OC1)
    has_epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts('C1OC1'))

    if not (has_hydroxyl or has_ketone or has_epoxide):
        return False, "No typical oxygen-containing functional groups found in xanthophyll"

    # If all checks passed, classify as xanthophyll
    return True, "Molecule is likely a xanthophyll (oxygenated carotenoid)"

# Example usage
# result, reason = is_xanthophyll("CC(\C=C\C=C(C)\C=C\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C[C@@H](O)CC1(C)C")
# print(result, reason)

# Note: Uncomment the example usage and replace the SMILES string with desired input to test the function.

# The function checks for the following:
# - Valid SMILES string
# - Sufficient number of carbons (typical carotenoids have around 40 carbons)
# - Presence of extensive conjugation (approximate estimation of conjugated double bonds)
# - Presence of oxygen atoms
# - Presence of typical functional groups in xanthophylls (hydroxyls, ketones, epoxides)

# If all criteria are met, the molecule is classified as a xanthophyll.
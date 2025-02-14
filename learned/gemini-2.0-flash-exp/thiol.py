"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a thiol group (-SH) attached to a carbon atom (aliphatic or aromatic)
    aliphatic_thiol_pattern = Chem.MolFromSmarts("[SH1][CX4]") # [SH1] explicitly captures sulfur with one H
    aromatic_thiol_pattern = Chem.MolFromSmarts("[SH1][c]")  # [SH1] explicitly captures sulfur with one H

    aliphatic_matches = mol.GetSubstructMatches(aliphatic_thiol_pattern)
    aromatic_matches = mol.GetSubstructMatches(aromatic_thiol_pattern)

    if not aliphatic_matches and not aromatic_matches:
        return False, "No thiol group (-SH) attached to a carbon atom found"

    return True, "Molecule contains a thiol group attached to a carbon"
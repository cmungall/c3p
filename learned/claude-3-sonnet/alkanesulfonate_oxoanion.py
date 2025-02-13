"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI:35518 alkanesulfonate oxoanion

An alkanesulfonate oxoanion is defined as an alkanesulfonate in which the carbon at position 1
is attached to R, which can represent hydrogens, a carbon chain, or other groups.
"""

from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alkanesulfonate oxoanion pattern ([*]-C-S(=O)(=O)[O-])
    alkanesulfonate_oxoanion_pattern = Chem.MolFromSmarts("[*]-C-S(=[O])(=[O])[O-]")
    if not mol.HasSubstructMatch(alkanesulfonate_oxoanion_pattern):
        return False, "Not an alkanesulfonate oxoanion"

    return True, "Molecule is an alkanesulfonate oxoanion"
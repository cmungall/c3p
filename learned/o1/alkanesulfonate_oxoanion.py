"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""

from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is an alkanesulfonate in which the carbon at position 1 is attached
    to R, which can represent hydrogens, a carbon chain, or other groups.
    It contains a sulfonate group (-S(=O)(=O)[O-]) attached to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alkanesulfonate oxoanion
    # Sulfonate group (-S(=O)(=O)[O-]) attached to a carbon
    pattern = Chem.MolFromSmarts('[C]-[S;D4;!R](=O)(=O)[O-]')

    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if the molecule contains the alkanesulfonate oxoanion pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains alkanesulfonate oxoanion group"
    else:
        return False, "Does not contain alkanesulfonate oxoanion group"
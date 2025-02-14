"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: aldehyde
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound with a carbonyl group bonded to at least one hydrogen atom and to one R group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldehyde functional group SMARTS pattern
    # Matches carbon with at least one hydrogen, double bonded to oxygen, single bonded to carbon or hydrogen
    aldehyde_pattern = Chem.MolFromSmarts("[$([CX3H][#6,#1])](=O)")

    # Check for aldehyde group
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an aldehyde functional group"
    else:
        return False, "Does not contain an aldehyde functional group"
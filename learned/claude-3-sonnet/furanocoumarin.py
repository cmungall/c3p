"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:32798 furanocoumarin

A furanocoumarin is any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways to give several isomers.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define furanocoumarin scaffold
    scaffold = Chem.MolFromSmarts("[o,O]1c2c(c3ccoc3)oc(=O)cc2c1")

    # Check if molecule contains furanocoumarin scaffold
    if mol.HasSubstructMatch(scaffold):
        return True, "Contains furanocoumarin scaffold (fused furan and coumarin rings)"
    else:
        return False, "Does not contain furanocoumarin scaffold"

    # Alternative approach using SMARTS pattern:
    # furanocoumarin_pattern = Chem.MolFromSmarts("[o,O]1c2c(c3ccoc3)oc(=O)cc2c1")
    # if mol.HasSubstructMatch(furanocoumarin_pattern):
    #     return True, "Contains furanocoumarin scaffold"
    # else:
    #     return False, "Does not contain furanocoumarin scaffold"
"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin. The fusion may occur in different ways,
    leading to several isomers (linear and angular).

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

    # Define SMARTS patterns for linear and angular furanocoumarins
    # Linear furanocoumarin (psoralen core)
    linear_smarts = 'O=C1C=CC2=CC=COC2=C1'
    linear_pattern = Chem.MolFromSmarts(linear_smarts)
    if linear_pattern is None:
        return False, "Invalid linear furanocoumarin SMARTS pattern"

    # Angular furanocoumarin (angelicin core)
    angular_smarts = 'O=C1C=CC2=COC=CC2=C1'
    angular_pattern = Chem.MolFromSmarts(angular_smarts)
    if angular_pattern is None:
        return False, "Invalid angular furanocoumarin SMARTS pattern"

    # Check for linear furanocoumarin core
    if mol.HasSubstructMatch(linear_pattern):
        return True, "Contains linear furanocoumarin core (psoralen type)"

    # Check for angular furanocoumarin core
    if mol.HasSubstructMatch(angular_pattern):
        return True, "Contains angular furanocoumarin core (angelicin type)"

    return False, "Does not contain furanocoumarin core"
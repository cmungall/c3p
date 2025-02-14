"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: Furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin structure.

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

    # Define the coumarin structure as a SMARTS pattern
    coumarin_pattern = Chem.MolFromSmarts("c1cc2ccc(=O)oc2c1")

    # Define the furan ring structure as a SMARTS pattern
    furan_pattern = Chem.MolFromSmarts("c1occc1")

    # Check presence of coumarin
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin structure found"
    
    # Check presence of furan
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"
    
    # Define a general pattern for furocoumarin fusion
    # This could be complex, involving various isomeric forms, but as a simple check:
    fusion_pattern = Chem.MolFromSmarts("c1occc1:c2ccccc2")
    if not mol.HasSubstructMatch(fusion_pattern):
        return False, "No furan-coumarin fusion ring found"

    return True, "Structure contains furan ring fused with a coumarin"
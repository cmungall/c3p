"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is defined as a nitrogen atom bonded to two hydrogen atoms and one other substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES, remove any salt info
    smiles = smiles.split(".")[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a primary amine (N with 2 hydrogens and 1 bond to any other atom)
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2]~[!H]")

    # Check for the substructure match
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains a primary amine"
    else:
        return False, "Does not contain a primary amine"
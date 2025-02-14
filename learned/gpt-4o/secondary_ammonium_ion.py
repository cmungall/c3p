"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has the general structure R2NH2+, formed by protonation of a secondary amine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the secondary ammonium ion pattern
    # N with two carbon attachments and a positive charge of 1
    secondary_ammonium_pattern = Chem.MolFromSmarts("[NH2+][C][C]")

    # Check if the molecule has the secondary ammonium ion structure
    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        return True, "Contains secondary ammonium ion pattern [NH2+][C][C]"
    else:
        return False, "Does not contain a secondary ammonium ion pattern"
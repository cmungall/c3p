"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion arises from the protonation of a secondary amine, indicated by
    a charged nitrogen ([+]) bound to exactly two distinct carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to match secondary ammonium ion: [NH2+] in context
    # Nitrogen carries a positive charge, two carbon neighbors
    secondary_ammonium_smart = Chem.MolFromSmarts("[NH2+][CH3C!H0][CH3C!H0]")

    # Identify if the molecule matches the secondary ammonium ion pattern.
    if mol.HasSubstructMatch(secondary_ammonium_smart):
        return True, "Contains protonated secondary amine group forming secondary ammonium ion"
    
    return False, "Does not contain the features of a secondary ammonium ion"
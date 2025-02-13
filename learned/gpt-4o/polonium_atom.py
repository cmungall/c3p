"""
Classifies: CHEBI:33313 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polonium atom, False otherwise
        str: Reason for classification
    """
    
    # Note: we typically expect SMILES of polonium atoms to be like "[210Po]" 
    # Check if the SMILES string is in an expected format for polonium
    if not smiles.startswith('[') or not smiles.endswith('Po]'):
        return False, "SMILES does not represent a polonium atom"

    # Try to extract the isotope number
    try:
        isotope_number = int(smiles[1:-3])
    except ValueError:
        return False, "Isotope number not recognized"
    
    # Polonium has isotopes in a specific range, though for this example we won't restrict
    # the exact values in this code.
    if isotope_number < 189 or isotope_number > 218:
        return False, f"Isotope number {isotope_number} not typical for known polonium isotopes"

    return True, f"SMILES represents polonium atom with isotope number {isotope_number}"
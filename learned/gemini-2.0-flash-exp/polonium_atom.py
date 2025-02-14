"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: CHEBI:33326 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.
    A polonium atom is represented by [nnnPo] where nnn is the atomic mass number.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polonium atom, False otherwise
        str: Reason for classification
    """
    if smiles is None:
        return False, "Input SMILES string is None"
    if not isinstance(smiles, str):
        return False, "Input is not a string"
    
    # Check for the correct format: [nnnPo]
    if not smiles.startswith('[') or not smiles.endswith(']'):
       return False, "Does not have the correct brackets"
    
    content = smiles[1:-1]
    if not content.endswith("Po"):
      return False, "Does not end with Po"
    
    mass_number = content[:-2]
    if not mass_number.isdigit():
       return False, "Mass number is not a digit"
    
    try:
        int(mass_number) # check if the mass number is actually a number
    except ValueError:
         return False, "Invalid mass number"
    
    return True, "The input is a polonium atom"
"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide contains a three-membered ring involving an oxygen and two carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Define an epoxide pattern: A three-membered ring with one oxygen atom
    epoxide_pattern = Chem.MolFromSmarts("[C]1OC1")
    
    # Check if the molecule contains an epoxide ring
    if mol.HasSubstructMatch(epoxide_pattern):
        return (True, "Contains a three-membered cyclic ether (epoxide)")

    return (False, "Does not contain a three-membered cyclic ether (epoxide)")
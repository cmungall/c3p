"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is a compound with a 1-benzopyran skeleton with an aryl substituent at position 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the flavonoid structural pattern
    # A flavonoid contains a 1-benzopyran ring system (C6-C3-C6)
    flavonoid_pattern = Chem.MolFromSmarts("c1cc2occc2c(c1)-c1ccccc1")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No 1-benzopyran structure with an aryl substituent at position 2 found"

    return True, "Contains 1-benzopyran structure with an aryl substituent at position 2"

# Example Usage:
# is_flavonoid("Oc1ccc(cc1O)-c1[o+]c2cc(O)c(O)c(O)c2cc1O")
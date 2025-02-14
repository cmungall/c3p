"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for steroid backbone (cyclopentanoperhydrophenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CC3C(C2)CCC4C3CCCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # SMARTS pattern for 3beta-hydroxy group
    hydroxy_3beta_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@H](CC2C(C1)CCC3C2CCC4C3CCCC4)')
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "3beta-hydroxy group not found"
    
    return True, "Contains steroid backbone and a 3beta-hydroxy group"
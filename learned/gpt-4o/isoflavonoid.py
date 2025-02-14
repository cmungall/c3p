"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for 1-benzopyran
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc1OCO") # Simplified pattern for core structure
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran backbone found"
        
    # Check for aryl group (i.e., ring) attached at position 3 of chromene
    aromatic_ring_pattern = Chem.MolFromSmarts("c1[cR2]oc2c1ccc(-c3ccccc3)c2=O") # Pattern for an aryl group at position 3
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aryl substituent at position 3 found"

    return True, "Contains 1-benzopyran with an aryl substituent at position 3"

# Testing a sample isoflavonoid SMILES
result, reason = is_isoflavonoid("COc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O")
print(result, reason)
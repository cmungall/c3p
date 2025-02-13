"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # A 1-benzopyran core: Benzene (c1ccccc1) fused to Pyran (O2c3c(C=CC=C3)CC(=C/1)\C2=O) with open valency at 4
    benzopyran_pattern = Chem.MolFromSmarts('O1C=CC=CC1C2=CC=CC=C2')
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "1-benzopyran core not found"
    
    # Check for aryl group attached at position 4 of the benzopyran core
    aryl_at_pos4_pattern = Chem.MolFromSmarts('O1C=CC=CC1-c2cccc3')
    if not mol.HasSubstructMatch(aryl_at_pos4_pattern):
        return False, "Aryl group not found at position 4"
    
    return True, "Contains a 1-benzopyran core with an aryl group at position 4"
"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:77834 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxy group at the 3-position in the alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]12CCC[C@H]3[C@@H]([C@@]1(CC[C@@H]2O)C)CCC4=CC(=O)CC[C@]34C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 3alpha-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@]([H])(O)([H])[C@@]12CCC[C@@H]3[C@@H]([C@@]1(CC[C@@H]2O)C)CCC(=O)C[C@]13C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3alpha-hydroxy group found"
    
    return True, "Contains a steroid backbone with a hydroxy group at the 3-position in the alpha configuration"
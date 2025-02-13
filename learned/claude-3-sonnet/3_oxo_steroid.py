"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:35693 3-oxo steroid
A 3-oxo steroid is any oxo steroid where an oxo substituent is located at position 3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]1(C[C@@]23[C@]([H])(CC[C@@]([H])(C2)C(C3)=O)CC[C@@]1([H])C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for oxo group at position 3
    oxo_pattern = Chem.MolFromSmarts("[C](=O)[C@]1(C)CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group at position 3"
    
    return True, "Contains a steroid backbone with an oxo group at position 3"
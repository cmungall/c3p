"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: CHEBI:38043 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is a steroid with a hydroxy group at position 11 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[#6,#8]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]" +
                                          "~[#6]~[#6]~[#6]~[#6]~[#6]~[#6,#8]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for 11-beta hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C@@]12[C@@H]([C@H]([C@@H]3[C@H]([C@@H]([C@H](C2)C)C4=CC(=O)CC[C@]34C)C)C)C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 11beta-hydroxy group found"
    
    return True, "Contains steroid backbone with 11beta-hydroxy group"
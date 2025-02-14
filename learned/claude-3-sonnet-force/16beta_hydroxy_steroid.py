"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:132026 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is a steroid with a hydroxy group at position 16 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]1(C[C@@H]2[C@@]1(CC[C@]1([C@]3([C@H](C[C@@H]4[C@@]3(CC[C@]2([C@@]4(C)C)C)C)C)C[C@@H](O)[C@@H]([C@@]1(C)C)O)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 16beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OH][C@]12[C@H](CC[C@@]1(C)C)CC[C@H]2[C@@H](O)C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 16beta-hydroxy group found"
    
    # Check molecular weight - steroids typically >200 Da
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for steroid"
    
    return True, "Contains steroid backbone with 16beta-hydroxy group"
"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: CHEBI:35595 11-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid has an oxo (=O) group at position 11 of the steroid core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@@H]3[C@@H]([C@@H]1CC[C@]2(C)O)CCC4=CC(=O)C=C45")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for oxo group at position 11
    oxo_pattern = Chem.MolFromSmarts("[C@@]1(CCC(=O)[C@@H]2[C@@]1(CCC3=CC(=O)C=C[C@@]23C)C)=O")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group at position 11"
    
    # Check molecular weight range (typical for steroids)
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, "Molecular weight out of typical steroid range"
    
    # Check for additional required features
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 19 or o_count < 2:
        return False, "Insufficient carbons or oxygens for steroid"
    
    return True, "Contains steroid backbone with oxo group at position 11"
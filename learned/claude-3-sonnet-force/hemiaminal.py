"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:27270 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic amino compound with an amino group and a hydroxy group
    attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for atoms with both amino (-NH2/-NH-) and hydroxy (-OH) groups
    hemiaminal_pattern = Chem.MolFromSmarts("[NH2,NH1][CH1][OH1]")
    if not mol.HasSubstructMatch(hemiaminal_pattern):
        return False, "No atoms with amino and hydroxy groups attached"
    
    # Check that amino and hydroxy groups are on same carbon
    hemiaminal_atoms = mol.GetSubstructMatches(hemiaminal_pattern)
    for match in hemiaminal_atoms:
        n_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        o_atom = mol.GetAtomWithIdx(match[2])
        if n_atom.GetNeighbors()[0].GetIdx() == c_atom.GetIdx() and o_atom.GetNeighbors()[0].GetIdx() == c_atom.GetIdx():
            return True, "Contains amino and hydroxy groups attached to the same carbon atom"
    
    return False, "Amino and hydroxy groups not attached to the same carbon atom"
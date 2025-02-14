"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:16676 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 1-benzopyran core
    benzopyran_pattern = Chem.MolFromSmarts("O=C1C=2C=CC=CC2=CC=C1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran core found"

    # Look for aryl substituent at position 3
    aryl_pattern = Chem.MolFromSmarts("[a]")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern, maxMatches=1, useChirality=True)
    if not aryl_matches:
        return False, "No aryl substituent found at position 3"
    
    aryl_idx = aryl_matches[0][0]
    atom = mol.GetAtomWithIdx(aryl_idx)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    if 2 not in neighbors:
        return False, "Aryl substituent not attached at position 3"

    return True, "Contains a 1-benzopyran core with an aryl substituent at position 3"
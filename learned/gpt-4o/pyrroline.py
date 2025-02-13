"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is defined as any organic heteromonocyclic compound with a structure based on a dihydropyrrole.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for different pyrroline patterns: Allow flexible positions for the double bond
    # Patterns cover 2-pyrroline, 3-pyrroline, etc.
    pyrroline_patterns = [
        Chem.MolFromSmarts("C1=CCNC1"),   # 2-pyrroline
        Chem.MolFromSmarts("C1=CNCC1"),   # 3-pyrroline
        Chem.MolFromSmarts("[n]=1[c,C][c,C][c,C][c,C]1"),  # Generic pyrroline with flexible atom type
    ]
    
    for pattern in pyrroline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a dihydropyrrole (pyrroline) ring"
    
    return False, "No pyrroline structure found"
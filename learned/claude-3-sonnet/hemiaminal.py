"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:35841 hemiaminal
A hemiaminal is any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.
"""

from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.

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
    
    # SMARTS pattern for hemiaminal substructure
    hemiaminal_pattern = Chem.MolFromSmarts("[NX3][CX4][OX2H]")
    
    # Check if the molecule contains the hemiaminal substructure
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains amino and hydroxy groups attached to the same carbon"
    else:
        return False, "No hemiaminal substructure found"
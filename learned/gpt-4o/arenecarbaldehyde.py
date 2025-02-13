"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde in which the carbonyl group is attached to an aromatic moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define aldehyde pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    
    # Check if the molecule contains an aldehyde group
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Find aldehyde matches
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    # Define aromatic ring pattern with attached aldehyde
    aromatic_aldehyde_pattern = Chem.MolFromSmarts("c-[CX3H1]=O")

    # Check if any match fits the aromatic pattern directly
    if mol.HasSubstructMatch(aromatic_aldehyde_pattern):
        return True, "Aromatic ring with attached aldehyde group found"
    
    return False, "No aromatic ring directly bound to an aldehyde group found"
"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:32879 arenecarbaldehyde
Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.

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

    # Detect aldehydes
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check for aromatic neighbors
    for match in aldehyde_matches:
        carbonyl_carbon_idx = match[0]
        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
        aromatic_neighbors = [
            neighbor for neighbor in carbonyl_carbon.GetNeighbors()
            if neighbor.GetIsAromatic()
        ]
        if aromatic_neighbors:
            return True, "Contains an aldehyde group attached to an aromatic ring"

    return False, "Aldehyde group not attached to an aromatic ring"
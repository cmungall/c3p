"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15677 primary alcohol
A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated 
carbon atom which has either three hydrogen atoms attached to it or only one other carbon atom
and two hydrogen atoms attached to it.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all hydroxy groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_atoms = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check each hydroxy to see if it is primary
    for hydroxy_idx in hydroxy_atoms:
        hydroxy_atom = mol.GetAtomWithIdx(hydroxy_idx)
        carbon = hydroxy_atom.GetNeighbors()[0]
        
        # Check if carbon is saturated (has 4 neighbors)
        if carbon.GetTotalDegree() != 4:
            continue
        
        # Check if carbon has 3 hydrogens OR 1 carbon and 2 hydrogens
        neighbors = [nbr.GetAtomicNum() for nbr in carbon.GetNeighbors()]
        if neighbors.count(1) == 3 or (neighbors.count(6) == 1 and neighbors.count(1) == 2):
            return True, "Contains a primary alcohol group (-OH attached to a saturated carbon with 3 H's or 1 C + 2 H's)"
    
    # If no primary alcohols found
    return False, "No primary alcohol groups found"
"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:33648 epoxide
An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule contains an epoxide group based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an epoxide group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for epoxide pattern ([O;r3])
    epoxide_pattern = Chem.MolFromSmarts("[O;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    
    if epoxide_matches:
        return True, "Contains a 3-membered cyclic ether (epoxide) group"
    else:
        return False, "No epoxide group found"
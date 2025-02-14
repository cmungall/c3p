"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:28722 isoflavone

An isoflavone is defined as any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavone(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isoflavone skeleton
    isoflavone_pattern = Chem.MolFromSmarts("c1cc(-c2coc3cc(O)ccc3c2=O)ccc1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Molecule does not contain isoflavone skeleton"
    
    # Check for substituents
    substituted_pattern = Chem.MolFromSmarts("c1ccc(c2coc3ccccc3c2=O)cc1")
    if mol.HasSubstructMatch(substituted_pattern):
        return True, "Molecule contains isoflavone skeleton with substituted derivatives"
    
    # Check for hydroxy and methoxy groups
    hydroxy_pattern = Chem.MolFromSmarts("O")
    methoxy_pattern = Chem.MolFromSmarts("OC")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    if hydroxy_matches or methoxy_matches:
        return True, "Molecule contains isoflavone skeleton with hydroxy or methoxy substituents"
    
    return True, "Molecule contains unsubstituted isoflavone skeleton"
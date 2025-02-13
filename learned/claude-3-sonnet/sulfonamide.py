"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:35526 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as an amide of a sulfonic acid: RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonamide pattern RS(=O)2NR'2
    sulfonamide_pattern = Chem.MolFromSmarts("[S+2]([N])([O-])(=[O])")
    matches = mol.GetSubstructMatches(sulfonamide_pattern)
    
    if matches:
        return True, "Contains sulfonamide group RS(=O)2NR'2"
    else:
        return False, "Does not contain sulfonamide group"
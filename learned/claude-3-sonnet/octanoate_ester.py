"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:36977 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule contains an octanoate (caprylic acid) ester group based on its SMILES string.
    An octanoate ester is any fatty acid ester in which the carboxylic acid component is octanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an octanoate ester group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for octanoate (caprylic acid) pattern (CCCCCCCC(=O)O)
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    
    # Check if pattern matches
    if mol.HasSubstructMatch(octanoate_pattern):
        return True, "Contains octanoate (caprylic acid) ester group"
    else:
        return False, "No octanoate ester group found"
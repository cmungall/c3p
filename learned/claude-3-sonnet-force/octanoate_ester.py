"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:35930 octanoate ester
Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for octanoate group (CCCCCCCC(=O)O-)
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)[O-]")
    octanoate_matches = mol.GetSubstructMatches(octanoate_pattern)
    
    # Look for ester bond (CCCCCCCC(=O)O-C)
    ester_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(octanoate_matches) > 0 or len(ester_matches) > 0:
        return True, "Contains octanoate group and ester bond"
    else:
        return False, "Does not contain octanoate group or ester bond"
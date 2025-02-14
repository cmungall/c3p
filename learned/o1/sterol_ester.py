"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:35366 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H0][#6]')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Look for steroid skeleton
    steroid_pattern_smiles = 'C1CC2CCC3CCCC4CCCCC4C3C2C1'  # Simplified steroid nucleus
    steroid_pattern = Chem.MolFromSmiles(steroid_pattern_smiles)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid skeleton found"
    
    # Check if ester group is connected to the steroid skeleton
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    
    # For each ester group, check if it's connected to the steroid skeleton
    for ester_match in ester_matches:
        ester_carbon_idx = ester_match[0]  # The carbonyl carbon
        ester_oxygen_idx = ester_match[1]  # The ester oxygen
        ester_carbon_attached_to_oxygen_idx = ester_match[2]  # The carbon attached to the ester oxygen
        
        # Check if the ester oxygen is connected to the steroid skeleton
        for steroid_match in steroid_matches:
            if ester_carbon_attached_to_oxygen_idx in steroid_match:
                return True, "Ester group connected to steroid skeleton"
    
    return False, "Ester group not connected to steroid skeleton"
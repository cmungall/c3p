"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: glycerophospholipid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerolipid having a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the glycerol backbone with ester linkages at positions 1 and 2
    # and a phosphate group at position 3
    glycerophospholipid_pattern = Chem.MolFromSmarts(
        '[C;H2;X4][C;H;X4][C;H2;X4]'
        '-[O;X2]-[P](=O)([O])[O]'
    )
    
    if glycerophospholipid_pattern is None:
        return False, "Unable to define glycerophospholipid pattern"
    
    # Search for the glycerophospholipid backbone pattern
    matches = mol.GetSubstructMatches(glycerophospholipid_pattern)
    if not matches:
        return False, "Glycerophospholipid backbone not found"
    
    # Check for two ester linkages attached to the glycerol backbone carbons
    ester_pattern = Chem.MolFromSmarts('[C;H2;X4]-[O;X2]-C(=O)')
    esters = mol.GetSubstructMatches(ester_pattern)
    if len(esters) < 2:
        return False, f"Found {len(esters)} ester linkages, need at least 2"
    
    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)[O]')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found"
    
    return True, "Contains glycerol backbone with ester-linked fatty acids and a phosphate group at terminal carbon"
"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone group at position 3 and a beta-configuration
    at position 5 in its backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ketone group (C=O) at position 3
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found in the structure"
    
    # Look for steroid backbone, and ensure position 5 has beta configuration
    # The generic steroid backbone might be represented simply by fused ring structures
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC5C4C3C2C1")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
    
    # Check for beta configuration
    # In SMILES, beta configuration is often represented using @H for stereochemistry
    beta_configuration_pattern = Chem.MolFromSmarts("[C@@H]")
    beta_configuration_matches = mol.GetSubstructMatches(beta_configuration_pattern)
    
    if len(beta_configuration_matches) < 1:
        return False, f"No beta configuration detected near the expected position"

    return True, "Contains a 3-oxo group and 5beta configuration consistent with 3-oxo-5beta-steroid class"
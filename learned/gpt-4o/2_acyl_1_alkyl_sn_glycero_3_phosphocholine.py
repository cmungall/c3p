"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    This class is characterized by a glycerol backbone with an alkyl ether at position 1, 
    an acyl ester at position 2, and a phosphocholine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    glycerol_backbone_pattern = Chem.MolFromSmarts("[C@H](CO)(COP(=O)(OCC[N+](C)(C)C)[O-])O")
    alkyl_ether_pattern = Chem.MolFromSmarts("COC")  # Simplified for alkyl ether presence check
    acyl_ester_pattern = Chem.MolFromSmarts("OC(=O)C")  # Generic acyl presence

    # Check for glycerol with phosphocholine at position 3
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "Structure doesn't match glycerol backbone with phosphocholine group"

    # Check for alkyl ether at position 1
    if not mol.HasSubstructMatch(alkyl_ether_pattern):
        return False, "No alkyl ether group found at position 1"

    # Check for acyl ester at position 2
    if not mol.HasSubstructMatch(acyl_ester_pattern):
        return False, "No acyl ester group found at position 2"

    return True, "Structure matches a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, with proper glycerol, alkyl, and acyl groups."
"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone with two hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"
    
    # Verify only one ester linkage (C(=O)OC) attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one ester linkage, found {len(ester_matches)}"
    
    # Verify presence of phosphate group directly linked to the glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not directly linked to glycerol backbone"
    
    # Check that there are no additional ester linkages indicating a diacylglycerol structure.
    # Specifically, ensure that no other acyl chains are attached
    additional_fatty_acid_pattern = Chem.MolFromSmarts("OC(=O)C")
    additional_fatty_acid_matches = mol.GetSubstructMatches(additional_fatty_acid_pattern)
    if len(additional_fatty_acid_matches) > 1:
        return False, "Additional acyl chain(s) found, indicating it may not be a lysophosphatidic acid"

    # All criteria met, it is a lysophosphatidic acid.
    return True, "The molecule is a lysophosphatidic acid"
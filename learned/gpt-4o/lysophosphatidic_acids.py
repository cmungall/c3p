"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid has a glycerol backbone, one acyl chain and a phosphate group.

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

    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for monoacyl group (-O-C(=O)-R)
    monoacyl_pattern = Chem.MolFromSmarts("OCOC(=O)")
    acyl_matches = mol.GetSubstructMatches(monoacyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"

    return True, "Contains glycerol backbone, one acyl group, and a phosphate group"
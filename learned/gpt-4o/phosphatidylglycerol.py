"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    Phosphatidylglycerol is characterized by a glycerol backbone, a phosphatidyl group, and fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible glycerol backbone pattern: capturing more variations
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)C")  # Simplified
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Broad phosphatidyl group recognizing flexibility of phosphatidic acid linkages
    phosphatidyl_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "Phosphatidyl group not properly detected"

    # Checking for two ester groups, indicative of fatty acid chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O")  # COOC representation of esters
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester linkages, found {len(ester_matches)}"
    
    return True, "Conforms to phosphatidylglycerol structure with glycerol backbone and phosphatidyl group"
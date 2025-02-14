"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol where one hydroxyl group is esterified with phosphoric acid and the other two with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone with two esterified hydroxyls and one phosphate ester
    glycerol_pattern = Chem.MolFromSmarts("OCC(O*)CO*)")
    phosphoric_acid_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    fatty_acid_ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Match glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Match phosphoric acid group
    if not mol.HasSubstructMatch(phosphoric_acid_pattern):
        return False, "No phosphoric acid group found"

    # Match two esterified fatty acid groups
    ester_matches = list(mol.GetSubstructMatches(fatty_acid_ester_pattern))
    if len(ester_matches) < 2:
        return False, "Less than 2 fatty acid ester groups found"

    # Check the chirality (if necessary in some cases, based on RS configuration)
    # This can be added as needed using RDKit stereochemistry features.

    return True, "Contains glycerol backbone with 2 fatty acid esters and a phosphoric acid ester"
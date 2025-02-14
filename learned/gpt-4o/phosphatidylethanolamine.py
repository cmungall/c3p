"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid with a phosphatidyl group esterified to ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with specific chiral center
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)C(O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with chiral center found"
    
    # Check for two ester groups (-O-C(=O)-) suggesting attachment of fatty acids
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for phosphoethanolamine group (P-OC-CN)
    phosphoethanolamine_pattern = Chem.MolFromSmarts("P(O)(OC)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    return True, "Contains glycerol backbone with 2 esterified fatty acids and a phosphoethanolamine group"
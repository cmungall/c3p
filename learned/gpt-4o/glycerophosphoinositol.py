"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with inositol esterified to the phosphate group at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True with reason if molecule is a glycerophosphoinositol, False otherwise with reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check glycerol backbone with O-P(=O)(O)O group at sn-3
    # The backbone must have a central carbon (CH) with two ester linkages, one leading to the phosphate (sn-3)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-phosphate backbone with correct linkage at sn-3 position found"
    
    # Ensure presence of inositol moiety attached to phosphate
    # This captures the inositol ring with primary and secondary hydroxyls
    inositol_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol moiety directly linked to the phosphate group found"
    
    # Check ester linkages indicating attachment of fatty acyl chains
    # Capture at least two ester bond patterns
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected at least 2 ester linkages for fatty acids, found {len(ester_matches)}"

    return True, "Contains glycerol-phosphate backbone with inositol moiety esterified at sn-3 position"
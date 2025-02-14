"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol phosphorylated at one or more hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a glycerol backbone. Example pattern: CCOC(=O)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)COC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-derived structure typical of phosphatidylinositol found"

    # Check for inositol ring: A six-membered ring with multiple hydroxy groups
    inositol_ring_pattern = Chem.MolFromSmarts("C1(O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_ring_pattern):
        return False, "No inositol ring found with hydroxy groups"

    # Check for phosphate groups attached to the inositol ring
    phosphate_group_pattern = Chem.MolFromSmarts("OP(=O)(O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1")
    phosphate_matches = mol.GetSubstructMatches(phosphate_group_pattern)

    if not phosphate_matches:
        return False, "No phosphate groups properly attached to the inositol"

    return True, "Contains glycerol backbone, inositol ring, and phosphorylated groups specific to phosphoinositides"
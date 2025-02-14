"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: CHEBI:26986 phosphatidylinositol phosphate

Any member of the phosphoinositide family of compounds, of which seven occur naturally.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate is a glycerolipid with a phosphorylated inositol head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Look for inositol ring pattern
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)OP(=O)([O-,O])[O-]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring with phosphate group found"

    # Count phosphate groups (could be 1-6)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-,O])[O-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    n_phosphates = len(phosphate_matches)
    if n_phosphates < 1 or n_phosphates > 6:
        return False, f"Found {n_phosphates} phosphate groups, expected 1-6"

    # Count fatty acid chains (could be different lengths)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    return True, "Contains glycerol backbone with 2 fatty acid chains and phosphorylated inositol head group"
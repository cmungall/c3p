"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: CHEBI:27526 glycerophosphoinositol
Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group at the sn-3 position of the glycerol backbone.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
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

    # Look for phosphate group (-O-P(=O)(-O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("O=P(-O)(-O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for inositol ring (C6H12O6 ring with 6 OH groups)
    inositol_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check that inositol is attached to phosphate at sn-3 position
    sn3_pattern = Chem.MolFromSmarts("CC(COP)CO")
    if not mol.HasSubstructMatch(sn3_pattern):
        return False, "Inositol not attached at sn-3 position"

    # Check for fatty acid chains (long carbon chains attached to other positions)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"

    return True, "Contains glycerol backbone with phosphoinositol attached at sn-3 and fatty acid chains at other positions"
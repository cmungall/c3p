"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has a myo-inositol ring, a glycerol backbone with two fatty acids, and a phosphate linking the inositol and glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a phosphatidylinositol, False otherwise
         str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for Inositol Ring
    inositol_pattern = Chem.MolFromSmarts("[C]1([C]([C]([C]([C]([C]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # 2. Check for Glycerol Backbone with phosphate and inositol.
    glycerol_phosphate_pattern = Chem.MolFromSmarts("C(COP(=O)(O)O[C]1([C]([C]([C]([C]([C]1O)O)O)O)O)O)") #glycerol with phosphate and inositol attached
    glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)

    if len(glycerol_phosphate_matches) == 0:
        glycerol_phosphate_pattern = Chem.MolFromSmarts("C(COP(=O)(O)[O][C]1([C]([C]([C]([C]([C]1O)O)O)O)O)O)") #glycerol with phosphate and inositol attached and an extra O
        glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if len(glycerol_phosphate_matches) == 0:
        return False, "No glycerol backbone with phosphate and inositol linkage found"

    # 3. Check for two fatty acid chains attached via esters.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if len(ester_matches) < 2:
      return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"


    # Check for number of rotatable bonds to exclude short chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8 :
        return False, "Chains too short to be fatty acids"

    # Count carbons and oxygens (approximate check).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if c_count < 20:
        return False, "Too few carbons for phosphatidylinositol"

    if o_count < 10:
      return False, "Too few oxygens for phosphatidylinositol"

    if p_count != 1:
        return False, f"Must have exactly 1 phosphorus, got {p_count}"

    return True, "Contains inositol, glycerol, phosphate, and two fatty acid chains"
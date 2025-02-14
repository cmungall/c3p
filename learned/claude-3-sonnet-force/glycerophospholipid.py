"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: CHEBI:18035 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerolipid with a phosphate group ester-linked to a
    terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
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

    # Look for phosphate group pattern (P(=O)(O)[O-,OH])
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-,OH])([O-,OH])")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for ester linkage between glycerol and phosphate
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Ensure the phosphate group is linked to a terminal carbon of glycerol
    for ester_match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(ester_match[1])
        if ester_atom.GetDegree() == 1:
            neighbor = ester_atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 8 and neighbor.IsInRingSize(3):
                continue  # Phosphate group is linked to glycerol backbone
            else:
                return False, "Phosphate group not linked to glycerol backbone"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    return True, "Contains glycerol backbone with phosphate group and fatty acid chains"
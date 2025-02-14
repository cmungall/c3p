"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: CHEBI:17754 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    Phosphatidylglycerols are glycerophospholipids where a phosphatidyl group is bound to glycerol.

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

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")  # Phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbones
    # Glycerol backbone: OCC(O)CO
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) < 2:
        return False, f"Found {len(glycerol_matches)} glycerol moiety(ies), need at least 2"

    # Check for phosphatidylglycerol headgroup
    # Phosphatidylglycerol headgroup: glycerol linked via phosphodiester bond to phosphate group
    pg_headgroup_pattern = Chem.MolFromSmarts("OCC(O)COP(O)(O)OCC(O)CO")
    if not mol.HasSubstructMatch(pg_headgroup_pattern):
        return False, "No phosphatidylglycerol headgroup found"

    # Check for ester-linked fatty acid chains
    # Ester linkage pattern: [C;H2,H1](OC(=O)C) - represents esterified glycerol carbons
    ester_linkage_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C](=O)[C]")
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Optional: Check if the fatty acid chains are of reasonable length (e.g., more than 8 carbons)
    # This step can be added if needed to improve specificity

    return True, "Molecule matches the phosphatidylglycerol structure"
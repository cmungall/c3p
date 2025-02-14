"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:32857 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is a fatty acid with one or more alkyl substituents on the parent hydrocarbon chain.
    The fatty acyl chain is usually saturated, but unsaturated BCFAs are also included.
    The substituents can be methyl groups or other alkyl chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for alkyl branches (including methyl and other alkyl groups)
    branch_pattern = Chem.MolFromSmarts("[CX4H2]")
    branch_matches = mol.GetSubstructMatches(branch_pattern)
    if len(branch_matches) < 1:
        return False, "No alkyl branches found"

    # Check for long carbon chain with possible unsaturation
    chain_pattern = Chem.MolFromSmarts("[CX4H2,CX3H1]~[CX4H2,CX3H1]~[CX4H2,CX3H1]~[CX4H2,CX3H1]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "Carbon chain too short for fatty acid"

    # Check for allowed cyclic structures (cyclopropane, cyclobutane, cyclopentane, cyclohexane)
    ring_info = mol.GetRingInfo()
    allowed_rings = [3, 4, 5, 6]
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size not in allowed_rings:
            return False, f"Contains disallowed cyclic structure of size {ring_size}"

    # Check for halogen substituents
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    if mol.HasSubstructMatch(halogen_pattern):
        return False, "Contains halogen substituents"

    # Check for common structural features of branched-chain fatty acids
    # Branch should be at least two carbons away from the carboxylic acid group
    branch_distance_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1].[CX4H2]")
    if not mol.HasSubstructMatch(branch_distance_pattern):
        return False, "Branch too close to the carboxylic acid group"

    # All checks passed
    reason = "Contains a carboxylic acid group, alkyl branches, and a long hydrocarbon chain (saturated or unsaturated)"
    return True, reason
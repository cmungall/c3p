"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoid monoterpenoids typically involve a cyclopentane ring fused to a six-membered oxygen heterocycle,
    or an open lactone form as seen in secoiridoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize major motif: cyclopentane ring
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")  # Basic cyclopentane
    if mol.HasSubstructMatch(cyclopentane_pattern):
        # Check for fusion with overlapping oxygen heterocycle
        oxygen_fusion_pattern = Chem.MolFromSmarts("[C;R]1[C;R][C;R][C;R]O1")
        if mol.HasSubstructMatch(oxygen_fusion_pattern):
            return True, "Found cyclopentane fused with oxygen-containing heterocycle"

    # Recognize if the structure is a possible secoiridoid (open form)
    results = mol.GetSubstructMatches(Chem.MolFromSmarts("O1CC=C[C;R][C;R][C;R]1"))
    if results:
        return True, "Found open lactone structure typical of secoiridoids"

    # Capture alternate overlapped rings scenarios (common in more complex iridoids)
    # e.g., 'dicyclohexyl' or 'methenyl-spiro' modifications, typical of secoiridoid variations
    alt_struct = Chem.MolFromSmarts("[C;R]1[C;R](O)[C;R][C;R](=O)[C;R]1")
    if mol.HasSubstructMatch(alt_struct):
        return True, "Typical iridoid structural motif (including secoiridoid content) recognized"

    return False, "Iridoid monoterpenoid characteristic structure not detected"
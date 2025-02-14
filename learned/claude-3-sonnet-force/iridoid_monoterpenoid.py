"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:35903 iridoid monoterpenoid

One of a class of monoterpenoids biosynthesized from isoprene and often intermediates in the biosynthesis of alkaloids. Iridoids usually consist of a cyclopentane ring fused to a six-membered oxygen heterocycle; cleavage of a bond in the cyclopentane ring gives rise to the subclass known as secoiridoids.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

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
    
    # Look for cyclopentane ring fused to a 6-membered oxygen heterocycle
    iridoid_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1]1[CR1]2[CR1][CR1][OR1][CR1][CR1]2")
    iridoid_match = mol.GetSubstructMatches(iridoid_pattern)
    if not iridoid_match:
        return False, "No iridoid core structure found"
    
    # Check if monoterpenoid (C10H16)
    formula = AllChem.CalcMolFormula(mol)
    if formula != "C10H16O":
        return False, "Not a monoterpenoid (formula does not match C10H16O)"
    
    # Check for secoiridoid subclass (cyclopentane ring cleaved)
    secoiridoid_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1]1[CR1]2[CR1][CR1][OR1][CR1][CR1]2[CR1]")
    secoiridoid_match = mol.GetSubstructMatches(secoiridoid_pattern)
    if secoiridoid_match:
        return True, "Iridoid monoterpenoid of the secoiridoid subclass"
    
    return True, "Iridoid monoterpenoid"
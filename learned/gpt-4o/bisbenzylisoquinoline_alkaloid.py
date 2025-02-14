"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determine if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
   
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the bisbenzylisoquinoline alkaloid pattern, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for two benzylisoquinoline units (simplistic heuristic pattern)
    benzylisoquinoline_pattern = Chem.MolFromSmarts('c1ccc2c(c1)CN3C=CC=C3C=C2')
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, "Does not contain two benzylisoquinoline units"

    # Check for ether bridges (-O-)
    ether_pattern = Chem.MolFromSmarts('[OD2]([#6])[#6]')
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) == 0:
        return False, "No ether bridges detected"

    return True, "Classified as bisbenzylisoquinoline alkaloid based on structural patterns"
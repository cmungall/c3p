"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    An organic sulfide has the structure R-S-R', where R and R' are organic groups (not hydrogen).
    The sulfur atom is bonded to two carbon atoms and no other atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the organic sulfide pattern: sulfur atom (S) with degree 2 bonded to two carbons (C)
    sulfide_pattern = Chem.MolFromSmarts('[S;D2]([#6])([#6])')
    if sulfide_pattern is None:
        return False, "Invalid sulfide SMARTS pattern"

    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(sulfide_pattern)
    if matches:
        return True, "Contains sulfur atom bonded to two carbon atoms (organic sulfide)"
    else:
        return False, "No sulfur atom bonded to two carbon atoms found"
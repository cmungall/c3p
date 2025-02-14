"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the anilide core SMARTS pattern
    anilide_pattern = Chem.MolFromSmarts("[c]1ccccc1[N]([CX3](=[OX1]))")

    # Search for the core pattern
    if not mol.HasSubstructMatch(anilide_pattern):
        return False, "No anilide core found"

    return True, "Anilide core detected"
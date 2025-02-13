"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: Nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a carbon atom triple-bonded to a nitrogen atom (C#N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitrile pattern: carbon triple bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    
    # Check for the nitrile pattern in the molecule
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile (C#N) group found"
    
    return True, "Contains a nitrile (C#N) group"
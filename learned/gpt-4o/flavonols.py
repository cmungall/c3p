"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core structure of flavonol (3-hydroxyflavone) in SMARTS notation
    flavonol_pattern = Chem.MolFromSmarts('Oc1cc(O)c2c(c1)oc(-c1ccccc1)c(=O)c2')
    
    # Check if molecule has the flavonol core structure
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol core structure found"
    
    # Additional check for common flavonol decorations (hydroxy groups)
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count hydroxyl groups, flavonols generally have multiple hydroxyl groups
    if len(hydroxy_matches) < 2:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxy_matches)}"

    return True, "Contains the core structure of a flavonol with sufficient hydroxyl groups"
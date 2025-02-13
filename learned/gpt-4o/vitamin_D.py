"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D is characterized as a seco-steroid with specific hydroxyl groups and often
    a conjugated triene system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the seco-steroid structure
    # Broadly, this could be characterized by the presence of a broken B-ring in a steroid nucleus
    seco_steroid_pattern = Chem.MolFromSmarts("C1C=C2CCCCC2=C1") # Simplified seco-steroid pattern
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No seco-steroid structure found"
    
    # Check for hydroxyl groups usually present at specific positions
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, typically require at least 2"

    # Check for conjugated triene system (flexible configuration based on vitamin D variants)
    conjugated_triene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C") # Approximation
    if not mol.HasSubstructMatch(conjugated_triene_pattern):
        return False, "No conjugated triene system found"

    # If all characteristics are matched, classify as vitamin D
    return True, "Molecule matches key structural features of vitamin D"
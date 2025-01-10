"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate with a specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible SMARTS pattern for a cyclohexane ring in inositols
    myo_inositol_pattern = Chem.MolFromSmarts("C1([OH])C([OH])C([OH])C([OH])C([OH])C1[OH]")
    
    # Check for the myo-inositol backbone
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol backbone found"
    
    # Check for the presence of any phosphate groups attached to the inositol
    phosphate_pattern = Chem.MolFromSmarts("[$(O[P](=O)(O)O)]") # Matches phosphate groups
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"
    
    # Get a count of phosphate groups attached to the myo-inositol backbone
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups attached to inositol backbone"

    return True, "Contains myo-inositol backbone with one or more phosphate groups"
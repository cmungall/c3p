"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine has a nitrogen atom connected to two non-heteroatom carbon atoms
    and not engaged in amide, imide, or other specific functional groups or rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for secondary amine (nitrogen attached to two carbons not involved in specific groups)
    sec_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([#6])[#6]")  # N with degree 3, one hydrogen, two carbon atoms
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")  # Match amide pattern where N is involved
    
    # Check for secondary amine and no amide group involvement
    if mol.HasSubstructMatch(sec_amine_pattern) and not mol.HasSubstructMatch(amide_pattern):
        return True, "Nitrogen is bonded to two carbon atoms and not part of an amide, indicative of a secondary amine"

    return False, "No nitrogen atom bonded to exactly two carbons in a secondary amine structure found"
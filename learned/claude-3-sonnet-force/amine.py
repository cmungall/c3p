"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:33877 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two or three
    hydrogen atoms by hydrocarbyl groups, or a nitrogen-containing heterocyclic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for different types of amines
    primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NH1]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NH0]")
    cyclic_amine_pattern = Chem.MolFromSmarts("[NX3]")
    
    # Check for matches
    primary_amine_matches = mol.GetSubstructMatches(primary_amine_pattern)
    secondary_amine_matches = mol.GetSubstructMatches(secondary_amine_pattern)
    tertiary_amine_matches = mol.GetSubstructMatches(tertiary_amine_pattern)
    cyclic_amine_matches = mol.GetSubstructMatches(cyclic_amine_pattern)
    
    if primary_amine_matches or secondary_amine_matches or tertiary_amine_matches or cyclic_amine_matches:
        reason = "Contains"
        if primary_amine_matches:
            reason += " primary amine(s)"
        if secondary_amine_matches:
            reason += " secondary amine(s)"
        if tertiary_amine_matches:
            reason += " tertiary amine(s)"
        if cyclic_amine_matches:
            reason += " cyclic amine(s)"
        return True, reason
    
    # Exclude quaternary nitrogen atoms
    quaternary_nitrogen_pattern = Chem.MolFromSmarts("[N+]")
    quaternary_nitrogen_matches = mol.GetSubstructMatches(quaternary_nitrogen_pattern)
    if quaternary_nitrogen_matches:
        return False, "Contains quaternary nitrogen atom(s)"
    
    # If no amine patterns found, return False
    return False, "No amine groups found"
"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are complex structures typically derived from triterpene skeletons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least 30 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Contains only {c_count} carbon atoms, requires at least 30 for typical triterpenoids"
    
    # Define skeletal structures typical of triterpenoids
    triterpenoid_patterns = [
        Chem.MolFromSmarts("C1CCC(C2)(C)CC(CC2C1)C"),  # Basic sterol skeleton
        Chem.MolFromSmarts("C1(CCCCC1)C2(CCCCC2)"),    # Pentacyclic triterpenoid-like
        # More patterns can be added to encompass diverse triterpenoid structures
    ]
    
    # Check against triterpenoid skeletal patterns
    for pattern in triterpenoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches triterpenoid-like skeletal pattern"

    # Check for common functional groups on a larger carbon skeleton
    common_functional_groups = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H]"),  # Carboxylic acid or ester
        Chem.MolFromSmarts("[OX2H][CX4]"),      # Alcohol
        Chem.MolFromSmarts("[NX3H1]")           # Amine, less common but possible
    ]
    
    if any(mol.HasSubstructMatch(fg) for fg in common_functional_groups):
        return True, "Contains functional groups reminiscent of triterpenoids"
    
    return False, "Does not match typical triterpenoid structures"
"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is classified as a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes, typically having a C10 skeleton that may be rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 16:
        return False, f"Carbon count of {c_count} falls outside typical monoterpenoid rearranged range (8-16)"

    # Check for common monoterpenoid functional groups, extended to include aldehyde and acid
    functional_groups = [
        Chem.MolFromSmarts("[CX3][OH]"),  # Alcohol
        Chem.MolFromSmarts("C(=O)[C]"),  # Ketone
        Chem.MolFromSmarts("C-O-C"),     # Ether
        Chem.MolFromSmarts("C(=O)O[C]"), # Ester
        Chem.MolFromSmarts("C(=O)O"),    # Acid
        Chem.MolFromSmarts("C=O")        # Aldehyde
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in functional_groups):
        return False, "No characteristic functional groups for monoterpenoids found (alcohol, ketone, ether, ester, acid, aldehyde)"

    # Check for typical monoterpenoid carbon backbone structure (acyclic or cyclic C10 core)
    typical_patterns = [
        Chem.MolFromSmarts("C1CCCCC1"),   # Cyclohexane core
        Chem.MolFromSmarts("C=CCCCC"),    # Acyclic C10 structure
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in typical_patterns):
        return False, "No typical monoterpenoid carbon skeleton found (cyclic/acyclic structures)"

    return True, "Contains a rearranged C10 skeleton or functional group characteristic of monoterpenoids."
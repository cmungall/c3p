"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoids with specific structural characteristics.

    Args:
        smiles (str): SMILES string of molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use more general patterns for flavonoid core that can account for variations
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("O=c1cc2ccccc2oc1"),  # Basic flavone pattern
        Chem.MolFromSmarts("Oc1cc2c(O)cc(O)cc2oc1"),  # More complex flavonoid core
    ]

    # Check for the flavonoid core
    if not any(mol.HasSubstructMatch(p) for p in flavonoid_core_patterns):
        return False, "Flavonoid core structure not found"

    # Check for hydroxyl groups - essential for anthoxanthins
    hydroxyl_group_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_group_pattern):
        return False, "No hydroxyl groups found, these are essential for anthoxanthins"

    # Check for glycosylation pattern (sugar components)
    sugar_pattern = Chem.MolFromSmarts("C1([OX2H][C@H]([O])C([O])[C@H](O)[C@H]1O)O")  # Simplified sugar pattern
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    has_sugar = len(sugar_matches) > 0
    
    # Return results based on identified patterns
    details = "Flavonoid core with hydroxyl groups"
    if has_sugar:
        details += " and glycosylation present"

    return True, details
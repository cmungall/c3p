"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26680 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols that occur in plants and have a core steroid structure, a hydroxyl at position 3 and side chains
    with some double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core steroid substructure (tetracyclic ring system) with a hydroxyl group at the 3 position, this is very important
    # Also the two methyl groups at positions 10 and 13 should be present
    steroid_pattern = Chem.MolFromSmarts("[C]1[C]([C])([H])[C]([H])([C]([H])([H])[C]2([C]([H])([H])[C]3([C]([H])([H])[C]([H])([C]([H])([H])[C]4([C]([H])([H])[C]([H])([C]([H])(O)[C]([H])([H])[C]14)[C]([H])([H])[H])[C]23)[H])([H])[H])")

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not have the core steroid structure"
    
    # Check for a hydroxyl group (-OH) on the core structure and confirm a single OH
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 1:
        return False, f"Must have exactly 1 hydroxyl group, found {len(hydroxyl_matches)}"

    # If all checks pass, it's likely a phytosterol
    return True, "Has core steroid structure with a single hydroxyl group"
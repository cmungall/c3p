"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.
    Quinic acid derivatives generally feature a cyclohexane ring with hydroxyl and ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid or relevant derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a cyclohexane ring
    cyclohexane_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane ring found"

    # Check for multiple hydroxy groups on the cyclohexane
    hydroxy_pattern = Chem.MolFromSmarts("[C;R1]O")
    num_hydroxy_matches = len(mol.GetSubstructMatches(hydroxy_pattern))
    if num_hydroxy_matches < 3:
        return False, f"Only {num_hydroxy_matches} hydroxy groups found, fewer than required for quinic acid"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Check for esterification, which could indicate derivatives
    ester_pattern = Chem.MolFromSmarts("OC=O[C;R1]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if any(ester_matches):
        return True, "Features consistent with esterified quinic acid"
    
    return True, "Core structure matches quinic acid"

# Example use case
smiles_example = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"
print(is_quinic_acid(smiles_example))
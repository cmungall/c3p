"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid has a hydroxyl group at the 11th position with beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for steroid core
    # Note: This is a simplified version and might not match all variations perfectly
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"

    # Define SMARTS pattern for 11-beta-hydroxy group
    # Considering carbon ring numbering, hydroxy group on 11th position with correct stereochemistry
    # Assuming convention: [C@@H](O) for beta configuration at 11th carbon
    beta_hydroxy_pattern = Chem.MolFromSmarts("[CH]1C[C@H](O)C(C)C1")
    
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "11-beta-hydroxy group not found or incorrect configuration"

    return True, "Contains 11-beta hydroxy group in correct configuration for steroids"

# Example usage
example_smiles = "[H][C@@]12CC[C@](O)(C(=O)CO)[C@]1(CO)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C"
print(is_11beta_hydroxy_steroid(example_smiles))
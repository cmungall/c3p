"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    A 2-oxo monocarboxylic acid has:
    - A carbonyl group at the alpha position.
    - A carboxylic acid functional group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for 2-oxo monocarboxylic acid: R-C(=O)-C(=O)O
    pattern = Chem.MolFromSmarts("C(=O)C(=O)O")
    # Check if the pattern is present in the molecule
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a 2-oxo monocarboxylic acid substructure"
    else:
        return False, "Does not contain a 2-oxo monocarboxylic acid substructure"

# Example usage
# smiles = "CC(=O)C(O)=O"  # pyruvic acid
# print(is_2_oxo_monocarboxylic_acid(smiles))
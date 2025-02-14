"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: CHEBI:26566 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine derivative based on its SMILES string.
    A 6-aminopurine derivative is any compound having 6-aminopurine (adenine) as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine derivative, False otherwise
        str: Reason for classification
    """

    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible SMARTS pattern for adenine (6-aminopurine)
    # This pattern searches for the purine ring with an amino group at position 6
    adenine_smarts = "c1nc(N)nc2ncnc12"
    adenine_mol = Chem.MolFromSmarts(adenine_smarts)
    if adenine_mol is None:
        return False, "Error in creating adenine substructure"

    # Perform substructure search
    if mol.HasSubstructMatch(adenine_mol):
        return True, "Molecule contains adenine (6-aminopurine) substructure"
    else:
        return False, "Molecule does not contain adenine (6-aminopurine) substructure"
"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is defined as an N-acyl-amino acid where the amino acid is glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted SMARTS pattern for N-acylglycine: acyl-N-C-C(=O)O
    n_acylglycine_pattern = Chem.MolFromSmarts("N([CX3]=[OX1])[CX4][CX3](=[OX1])[OX2H1]")  # Key pattern for glycine-derived N-acyl compound
    
    # Check for the substructure match
    if mol.HasSubstructMatch(n_acylglycine_pattern):
        return True, "Contains N-acylglycine structure"

    return False, "Does not match N-acylglycine substructure"
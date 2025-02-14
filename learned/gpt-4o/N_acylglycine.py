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

    # Updated SMARTS pattern for N-acylglycine: acyl-NH-CH2-CO-O
    n_acylglycine_pattern = Chem.MolFromSmarts("N[C;H2][CX3](=O)[O;H1]") 
    if not mol.HasSubstructMatch(n_acylglycine_pattern):
        return False, "Does not match N-acylglycine substructure"

    # To ensure it's N-acyl, not just simple NH linkage
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")  
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Does not have N-acyl linkage"

    return True, "Contains N-acylglycine structure"
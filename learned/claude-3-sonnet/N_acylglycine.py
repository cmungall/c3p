"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI:48215 N-acylglycine
An N-acyl-amino acid in which amino acid specified is glycine.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine has an acyl group (R-C(=O)-) attached to the nitrogen of glycine (NH-CH2-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the N-acylglycine substructure pattern
    n_acylglycine_pattern = Chem.MolFromSmarts("[C;$(C(=O)N[CH2;$(CC(=O)O)])]")
    
    # Check if the molecule matches the pattern
    match = mol.GetSubstructMatch(n_acylglycine_pattern)
    if not match:
        return False, "No N-acylglycine substructure found"
    
    return True, "Contains an acyl group (R-C(=O)-) attached to the nitrogen of glycine (NH-CH2-COOH)"
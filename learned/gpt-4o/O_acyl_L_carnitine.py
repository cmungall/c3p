"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is an O-acylcarnitine in which the carnitine component has L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern for L-carnitine esters: Check for specific acyl linkage and L-carnitine stereochemistry
    acyl_carnitine_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C")
    # Inverse pattern for cases with 'C@@H' which is L-carnitine
    inverse_acyl_carnitine_pattern = Chem.MolFromSmarts("C(=O)O[C@@H](CC(=O)[O-])C[N+](C)(C)C")
    
    # Check for presence of either configuration
    if mol.HasSubstructMatch(acyl_carnitine_pattern):
        return True, "Correctly matches O-acyl-L-carnitine pattern"
    if mol.HasSubstructMatch(inverse_acyl_carnitine_pattern):
        return True, "Correctly matches O-acyl-L-carnitine pattern with inverse stereochemistry"
    
    return False, "Does not match O-acyl-L-carnitine pattern"
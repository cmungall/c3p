"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is characterized by a -C(=O)S- linkage to a coenzyme A structure.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found (-C(=O)S-)"
    
    # Define a SMARTS pattern for the coenzyme A backbone
    # This is a simplification capturing the phosphopantetheine with typical CoA adjacencies
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](O[C@H]1OP(=O)(O)O)ncnc2c(N)ncnc12")

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found (missing typical CoA structure)"
    
    # Check for additional structural requirements or other logical checks as necessary (e.g., molecular weight)

    return True, "Molecule contains thioester linkage to a coenzyme A moiety"
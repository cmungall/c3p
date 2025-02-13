"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is defined by having a phosphatidyl group esterified to 
    one of the hydroxy groups of inositol, a hexahydroxy cyclohexane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the inositol ring pattern: a six-carbon ring with all hydroxyl groups
    inositol_pattern = Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C(O)O1)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"
    
    # Identify the phosphatidyl group pattern: phosphate attached to glycerol
    phosphatidyl_pattern = Chem.MolFromSmarts("C(COP(=O)(O)O)CO")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl group found"

    # Check if the two patterns are connected correctly
    connected_pattern = Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C(O)O1)O-[P](=O)(O)-O-[C]-[C]-[O]-[C]-[O]")
    if not mol.HasSubstructMatch(connected_pattern):
        return False, "Inositol and phosphatidyl groups not connected correctly"
    
    return True, "The molecule is a phosphatidylinositol with correct linkages"
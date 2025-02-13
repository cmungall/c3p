"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:35914 acetate ester
An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for acetate group pattern (-O-C(=O)C)
    acetate_pattern = Chem.MolFromSmarts("[OX2]C(=O)C")
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)
    
    # Look for ester pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check if any ester group is also an acetate group
    for ester_match in ester_matches:
        if any(ester_match in acetate_match for acetate_match in acetate_matches):
            return True, "Contains an acetate ester group (-O-C(=O)C)"

    return False, "No acetate ester group found"
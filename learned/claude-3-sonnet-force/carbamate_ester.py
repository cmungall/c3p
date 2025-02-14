"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:35703 carbamate ester
Any ester of carbamic acid or its N-substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbamate functional group pattern (-O-C(=O)-N-)
    carbamate_pattern = Chem.MolFromSmarts("[OX2]C(=O)N")
    carbamate_matches = mol.GetSubstructMatches(carbamate_pattern)
    
    if len(carbamate_matches) > 0:
        return True, "Contains carbamate functional group (-O-C(=O)-N-)"
    else:
        return False, "Does not contain carbamate functional group"
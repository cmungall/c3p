"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate contains the functional group N=C=S.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for isothiocyanate group pattern N=C=S
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        return True, "Contains isothiocyanate group (N=C=S)"
    else:
        return False, "No isothiocyanate group (N=C=S) found"

# Example usage:
# result, reason = is_isothiocyanate("N(=C=S)CCCCS(=O)C")
# print(result, reason)
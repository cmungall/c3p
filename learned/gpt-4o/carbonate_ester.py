"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the general structure O=C(OX)OX, where X is any organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carbonate ester pattern: O=C(O[C])O[C]
    carbonate_ester_pattern = Chem.MolFromSmarts("O=C(O[*])O[*]")
    if not mol.HasSubstructMatch(carbonate_ester_pattern):
        return False, "No carbonate ester pattern found"

    return True, "Contains carbonate ester structure"

# Example usage
# print(is_carbonate_ester("COC(=O)OC"))  # Example for dimethyl carbonate
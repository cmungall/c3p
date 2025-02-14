"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha-hydroxy ketone pattern: [OH] next to the [C=O] connected by [C] with [H][R]
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts('[CX4H1][CX3](=O)[CX4](O)')
    if not mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return False, "No secondary alpha-hydroxy ketone structure found"

    return True, "Contains secondary alpha-hydroxy ketone structure"

# Example usage: Test the function with a known secondary alpha-hydroxy ketone.
smiles_example = "OC[C@@H]1CC=C(O)C(=O)[C@H]1O" # 5D-(5/6)-2,6-dihydroxy-5-(hydroxymethyl)cyclohex-2-en-1-one
result, reason = is_secondary_alpha_hydroxy_ketone(smiles_example)
print(f"Result: {result}, Reason: {reason}")
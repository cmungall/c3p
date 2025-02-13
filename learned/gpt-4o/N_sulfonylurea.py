"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas are characterized by a sulfonyl group attached to a nitrogen of a urea moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated pattern for N-sulfonylurea moiety
    # More flexible arrangement to account for common variations
    # N-sulfonylureas have a general form comprising an N-sulfonyl group (N-S(=O)(=O)) attached to a urea (N-C(=O)-N)
    nsulfonylurea_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3][S](=O)(=O)[#6]")
    alternative_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[SX4](=O)(=O)[NX3]")  # Consider reversal in structures

    # Check if the molecule has these patterns
    match1 = mol.HasSubstructMatch(nsulfonylurea_pattern)
    match2 = mol.HasSubstructMatch(alternative_pattern)

    if match1 or match2:
        return True, "Contains N-sulfonylurea moiety"
    else:
        return False, "Does not contain N-sulfonylurea moiety"
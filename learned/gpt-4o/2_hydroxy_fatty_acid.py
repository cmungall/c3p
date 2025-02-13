"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy functional group in the alpha- or 2-position
    relative to the carboxylic acid group, with a sufficiently long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for 2-hydroxy fatty acids
    # Match a carboxylic acid group (C(=O)O) with an adjacent hydroxy group (O)
    # This pattern includes both anionic and neutral forms
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[C;H1,H2]([O])C(=O)[O;H1,H0]")
    
    # Check for the correct substructure in the molecule
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No hydroxy group in the alpha position to the carboxylic acid group"
    
    # Ensuring the presence of a sufficient length aliphatic chain
    # Check for at least 8 consecutive carbon atoms indicating a longer chain
    chain_length_pattern = Chem.MolFromSmarts("[C][C][C][C][C][C][C][C]")
    if not mol.HasSubstructMatch(chain_length_pattern):
        return False, "Insufficient carbon chain length typical for a fatty acid"
    
    return True, "Contains hydroxy group at 2-position with appropriate aliphatic chain and carboxylic acid group"

# Testing examples
# Example that should match
print(is_2_hydroxy_fatty_acid("CCCCCCCC(O)C(=O)O"))  # Should return True
# Example that should not match
print(is_2_hydroxy_fatty_acid("OC(C(=O)O)C"))  # Should return False due to chain length
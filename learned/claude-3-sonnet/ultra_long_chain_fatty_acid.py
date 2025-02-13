"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:51165 ultra-long-chain fatty acid
An ultra-long-chain fatty acid is any very long-chain fatty acid which has a chain length greater than C27.
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Does not contain a carboxylic acid group"
    
    # Find the length of the longest carbon chain
    chain_length = 0
    for length in range(28, 100):  # Iterate over chain lengths from 28 to 99
        chain_pattern = Chem.MolFromSmarts(f"[C]~[C]~[C]~[C]~[C]~[C]~[C]~{f'~[C]' * (length - 7)}")
        if mol.HasSubstructMatch(chain_pattern):
            chain_length = length
            break
    
    if chain_length <= 27:
        return False, f"Carbon chain length is only {chain_length}, need greater than 27"
    
    # Additional checks (optional)
    # Check for the absence of rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structures, expected an acyclic compound"
    
    # Check for the presence of only single bonds in the carbon chain
    # ...
    
    # Check for the absence of additional functional groups
    # ...
    
    return True, f"Contains a carbon chain of length {chain_length} with a carboxylic acid group"
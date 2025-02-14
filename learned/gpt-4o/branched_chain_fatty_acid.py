"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A BCFA is characterized by the presence of a carboxylic acid group and one or more branches
    in its hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a BCFA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Define branching pattern as a methyl group attachment
    branch_pattern = Chem.MolFromSmarts("[C;D3]([C;D1])")  # Tertiary carbon with a methyl group
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No suitable branching pattern found"
    
    # Ensure branching patterns are not part of complex ring systems, or large non-linear structures
    for match in mol.GetSubstructMatches(branch_pattern):
        atom = mol.GetAtomWithIdx(match[0])
        if not atom.IsInRing():
            return True, "Contains a branched hydrocarbon chain with carboxylic acid group"
    
    return False, "Branches found are part of non-fatty systems or rings"

# Testing examples
smiles_examples = [
    "OC(=O)C=C(C)C",  # Expected: True (branched-chain fatty acid)
    "CCCCCCCCCCCCCCCCCC(=O)O",  # Expected: False (straight-chain fatty acid)
]

# Output results for testing examples
results = {sm: is_branched_chain_fatty_acid(sm) for sm in smiles_examples}
for sm, result in results.items():
    print(f"SMILES: {sm} -> {result}")
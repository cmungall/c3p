"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:27283 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count number of carboxylic acid groups - should be exactly one
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    if carboxyl_matches > 1:
        return False, "Multiple carboxylic acid groups found"
        
    # Check for aromatic rings - should have none
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
        
    # Count carbons (excluding the carboxyl carbon)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    chain_carbons = carbon_count - 1  # Subtract the carboxyl carbon
    
    if chain_carbons >= 6:
        return False, f"Chain too long ({chain_carbons} carbons excluding carboxyl group)"
    if chain_carbons < 1:
        return False, "No carbon chain found"
        
    # Check for non-hydrocarbon atoms (except the carboxylic acid oxygens)
    allowed_atoms = {6, 1, 8}  # C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains non-hydrocarbon atoms other than carboxylic acid"
            
    # Count oxygens - should be exactly 2 (from carboxylic acid)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count > 2:
        return False, "Contains oxygen atoms other than carboxylic acid"
        
    # All checks passed
    return True, f"Aliphatic monocarboxylic acid with {chain_carbons} carbons in chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27283',
        'name': 'short-chain fatty acid',
        'definition': 'An aliphatic monocarboxylic acid with a chain length of less than C6. '
                     'If any non-hydrocarbon substituent is present, the compound is not normally '
                     'regarded as a short-chain fatty acid.',
    }
}
"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.

    A long-chain fatty acid is defined as a fatty acid with a carbon chain length from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count total number of carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Ensure carbon chain length is between C13 and C22
    if carbon_count < 14 or carbon_count > 23:
        return False, f"Total carbon count {carbon_count} is not within the long-chain fatty acid range (14-23 considering COOH group)"

    # Check for a ring structure which is typically not part of long-chain fatty acids
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring(s), not a typical structure for long-chain fatty acids"

    return True, f"Valid long-chain fatty acid with {carbon_count} carbon atoms"
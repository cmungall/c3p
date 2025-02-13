"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid with a chain length of less than C6,
    with no non-hydrocarbon substituents outside the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carboxylic acid group: -C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count >= 6:
        return False, f"Contains {c_count} carbons, expected fewer than 6"
    
    # Check for non-hydrocarbon substituents (excluding carboxylic group itself)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {'H', 'C', 'O'}:
            return False, f"Contains non-hydrocarbon substituent: {atom.GetSymbol()}"
    
    # Ensure there are no aromatic atoms, i.e., aliphatic
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic structures, expected aliphatic"
    
    return True, "Valid short-chain fatty acid with carboxylic acid group and <6 carbon atoms"

# Example call 
# This should return (True, "Valid short-chain fatty acid with carboxylic acid group and <6 carbon atoms")
# print(is_short_chain_fatty_acid("CCCC(O)=O"))  # Represents butyric acid
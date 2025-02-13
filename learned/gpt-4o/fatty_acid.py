"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
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

    # Count carbon atoms to evaluate aliphatic chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Natural fatty acids commonly range from C4 to C28 in aliphatic count
    if not (4 <= carbon_count <= 28):
        return False, f"Carbon chain length of {carbon_count}, must be between 4 and 28"

    # Check for the presence of aromaticity
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings, not a fatty acid"

    # Check if the molecule is predominantly linear
    cyclic = any(atom.IsInRing() for atom in mol.GetAtoms())
    if cyclic:
        # Some natural fatty acids could have minimal cyclic systems (rare)
        return False, "Has cyclic structures which are typically not part of standard fatty acids"

    # Verify functionalities like hydroxyls or epoxides are accommodated
    # These rules can be expanded based on detailed chemical understanding
    
    return True, "Contains aliphatic carbon chain with a carboxylic acid group, fits fatty acid profile"
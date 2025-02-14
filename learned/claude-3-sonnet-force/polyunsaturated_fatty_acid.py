"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid contains more than one double bond and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, "Less than two double bonds found"
    
    # Check if double bonds are separated by at least one methylene group
    chain = Chem.MolToSmiles(mol)
    double_bond_indices = [match[0] for match in double_bond_matches]
    for i in range(len(double_bond_indices) - 1):
        start = double_bond_indices[i] + 2  # Skip double bond atoms
        end = double_bond_indices[i+1]
        if "C=C" in chain[start:end]:
            return False, "Double bonds are not separated by methylene groups"
    
    return True, "Contains more than one double bond and a carboxylic acid group"
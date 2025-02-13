"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: CHEBI:36975 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, "Less than two double bonds found"
    
    # Check for long carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2X4,CH3X3]~[CH2X4,CH3X3]~[CH2X4,CH3X3]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(carbon_chain_matches) < 3:
        return False, "Carbon chain too short for a fatty acid"
    
    return True, "Molecule contains a carboxylic acid group, more than one double bond, and a long carbon chain"
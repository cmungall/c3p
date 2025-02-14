"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:36195 trienoic fatty acid
A polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Count double bonds
    num_double_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_double_bonds != 3:
        return False, f"Found {num_double_bonds} double bonds, expected 3"

    # Check for long carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "No long carbon chain found"

    # Check for fatty acid
    fatty_acid_pattern = Chem.MolFromSmarts("CCC(=O)O")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Not a fatty acid"

    return True, "Contains a long carbon chain with three double bonds and a carboxylic acid group"
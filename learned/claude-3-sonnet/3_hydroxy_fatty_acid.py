"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:38120 3-hydroxy fatty acid

A 3-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the 
beta- or 3-position. beta-Hydroxy fatty acids accumulate during cardiac hypoxia,
and can also be used as chemical markers of bacterial endotoxins.
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid functional group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for hydroxy group at 3-position
    hydroxy_pattern = Chem.MolFromSmarts("[C;H3]([OH])([C;H2])[C;H2]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group at 3-position found"

    # Check for long aliphatic chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "No long aliphatic chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = AllChem.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    # All criteria satisfied, it's a 3-hydroxy fatty acid
    return True, "Contains a carboxylic acid group, a hydroxy group at the 3-position, and a long aliphatic chain"
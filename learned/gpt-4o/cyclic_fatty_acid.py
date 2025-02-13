"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is any fatty acid containing anywhere in its structure a ring of atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found, not a fatty acid"

    # Check for long carbon chain (minimum criteria could be, e.g., 6 continuous carbons)
    chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient carbon chain found, not a fatty acid"

    # Detect any cyclic structure
    if not mol.GetRingInfo().NumRings():
        return False, "No cyclic structure found"

    return True, "Contains cyclic structure and fatty acid characteristics"
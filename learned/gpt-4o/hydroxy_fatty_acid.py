"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is defined as a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("O=C[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for at least one hydroxy group (-O[H])
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy (OH) group found"
    
    # Check the carbon chain length is sufficient
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:  # Assume a practical minimum length for a hydroxy fatty acid
        return False, "Too few carbon atoms for a fatty acid"
    
    # Fatty acids are primarily linear
    # We allow some degree of non-linearity (branching/rings) but should be limited
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, "Too many rings detected, inconsistent with typical fatty acid structure"

    return True, "Molecule contains a carboxylic acid group and one or more hydroxy groups; matches hydroxy fatty acid structure"
"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:33594 tricarboxylic acid
An oxoacid containing three carboxy groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group pattern (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Check if there are exactly 3 carboxyl groups
    if len(carboxyl_matches) != 3:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, need exactly 3"

    # Check for oxoacid pattern (-C(=O)-OH)
    oxoacid_pattern = Chem.MolFromSmarts("C(=O)O")
    oxoacid_match = mol.HasSubstructMatch(oxoacid_pattern)
    if not oxoacid_match:
        return False, "Not an oxoacid (missing -C(=O)-OH group)"

    return True, "Contains three carboxyl groups (-COOH)"
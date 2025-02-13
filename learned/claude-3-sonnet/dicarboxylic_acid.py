"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35652 dicarboxylic acid
Any carboxylic acid containing two carboxy groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains two carboxyl groups (-C(=O)O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group pattern (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    # Check if there are exactly 2 carboxyl groups
    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, need exactly 2"

    # Check for connectivity between carboxyl groups
    # (i.e., they are not part of separate fragments)
    fragment_info = AllChem.GetMolFrags(mol)
    if len(fragment_info) > 1:
        return False, "Carboxyl groups are not connected"

    return True, "Contains two connected carboxyl groups (-C(=O)O)"
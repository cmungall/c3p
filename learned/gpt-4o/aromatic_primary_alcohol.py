"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl (-OH) group attached to a primary carbon that is
    directly bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS query for aromatic primary alcohol
    # Primary carbon bonded to -OH and directly connected to aromatic ring
    aromatic_primary_alcohol_smarts = Chem.MolFromSmarts('[c]-[CH2]-[OH]')

    # Check if the molecule matches the aromatic primary alcohol SMARTS pattern
    if mol.HasSubstructMatch(aromatic_primary_alcohol_smarts):
        return True, "Matches aromatic primary alcohol pattern"

    return False, "Does not match aromatic primary alcohol pattern"
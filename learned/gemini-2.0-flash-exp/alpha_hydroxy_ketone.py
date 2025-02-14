"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone with a hydroxy group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Define SMARTS pattern for a ketone
    ketone_pattern = Chem.MolFromSmarts("[CX3]=O")
    
    # 3. Find all ketone matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # 4. If there is no ketone, return False.
    if not ketone_matches:
      return False, "No ketone found."

    # 5. Define SMARTS pattern for alpha-hydroxy ketone
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX4](-[OX2H1])-[CX3](=O)")

    # 6. Check if an alpha-hydroxy ketone is found
    matches = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)

    if matches:
       return True, "Hydroxy group found on an alpha carbon of a ketone"
    else:
       return False, "No hydroxy group found on any alpha carbon of a ketone"
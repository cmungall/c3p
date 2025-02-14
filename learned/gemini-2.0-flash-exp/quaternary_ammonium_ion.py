"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is a positively charged nitrogen with four single bonds to non-hydrogen atoms (usually carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a quaternary ammonium ion
    # [N+] is a positively charged nitrogen
    # [C] means a carbon atom
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+]([C])([C])([C])[C]")
    
    if mol.HasSubstructMatch(quaternary_ammonium_pattern):
         return True, "Quaternary ammonium ion pattern found: Positively charged nitrogen with four single bonds to carbon atoms."
    
    return False, "No quaternary ammonium ion pattern found: Positively charged nitrogen with four single bonds to carbon atoms."
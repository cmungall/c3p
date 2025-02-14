"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:35671 chalcones
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is an aromatic ketone that consists of two aromatic rings
    joined by a three-carbon α,β-unsaturated carbonyl system (Ar-CH=CH-C(=O)-Ar).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the chalcone core SMARTS pattern
    # Pattern represents:
    # Aromatic ring - C=C - C=O - Aromatic ring
    # where the C=C is an alpha,beta-unsaturated ketone
    chalcone_pattern = Chem.MolFromSmarts("[a]-[#6]=[#6]-C(=O)-[a]")
    if chalcone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for the chalcone substructure
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Contains chalcone core structure"
    else:
        return False, "Does not contain chalcone core structure"
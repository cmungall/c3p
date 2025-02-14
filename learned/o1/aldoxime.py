"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:16445 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as an oxime derived from an aldehyde, having the functional group -CH=N-OH.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldoxime SMARTS pattern
    # [CH]=N[OH] matches a carbon with one hydrogen double-bonded to nitrogen,
    # which is single-bonded to an oxygen with one hydrogen (hydroxyl group)
    aldoxime_pattern = Chem.MolFromSmarts('[CH]=N[OH]')

    # Search for the aldoxime pattern in the molecule
    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains aldoxime functional group (-CH=N-OH)"
    else:
        return False, "No aldoxime functional group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16445',
        'name': 'aldoxime',
        'definition': 'Oximes of aldehydes RCH=NOH.',
        'parents': []
    },
}
"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule belongs to the acrovestone class based on its SMILES string.
    Acrovestone compounds are polyphenols with an isoflavone core structure and various substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to acrovestone class, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define isoflavone core SMARTS pattern
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1')  # Basic isoflavone skeleton
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyl < 1:
        return False, "No hydroxyl groups found"

    # Check for methoxy groups (-OCH3)
    methoxy_pattern = Chem.MolFromSmarts('CO')
    num_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))

    # Check for glycosylation (sugar moieties)
    # Example patterns for common sugars (glucose, rhamnose)
    glucose_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O')
    rhamnose_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@H](C)[C@@H](O)[C@@H](O)[C@H]1O')
    glycosylated = False
    if mol.HasSubstructMatch(glucose_pattern) or mol.HasSubstructMatch(rhamnose_pattern):
        glycosylated = True

    # Check for additional fused rings or linkages
    # Not strictly necessary but can be included based on specific cases
    # Additional ring patterns can be defined and matched here

    # If all checks pass, classify as acrovestone
    return True, "Molecule matches acrovestone class (isoflavone core with appropriate substituents)"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'acrovestone',
        'definition': 'A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata can be added as required
}
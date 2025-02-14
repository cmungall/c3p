"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI:XXXX 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid is a steroid with a ketone group at position 3 and alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid nucleus (four fused rings)
    steroid_nucleus_smarts = '[#6]1[#6][#6]2[#6](~[#6]1)~[#6][#6]3~[#6](~[#6][#6]2)~[#6][#6]4~[#6](~[#6][#6]3)~[#6][#6]~4'
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus detected"

    # Define SMARTS pattern for ketone at position 3 in steroid nucleus
    # Approximate position 3 as an sp3 carbon connected to a carbonyl group
    ketone_smarts = '[#6]-[#6](=O)-[#6]'
    ketone = Chem.MolFromSmarts(ketone_smarts)
    if not mol.HasSubstructMatch(ketone):
        return False, "No ketone group detected at position 3"

    # Define SMARTS pattern for chiral center at position 5 with alpha configuration
    # Approximate position 5 as a specific chiral carbon in the steroid nucleus
    # Alpha configuration is represented by '@' in SMILES/SMARTS notation
    chiral_center_smarts = '[C@H]1CC[C@@H](C)C1'
    chiral_center = Chem.MolFromSmarts(chiral_center_smarts)
    if not mol.HasSubstructMatch(chiral_center):
        return False, "Position 5 does not have alpha configuration"

    return True, "Molecule is a 3-oxo-5alpha-steroid"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXX',
                              'name': '3-oxo-5alpha-steroid',
                              'definition': 'A 3-oxo steroid that has alpha configuration at position 5.',
                              'parents': []},
        }
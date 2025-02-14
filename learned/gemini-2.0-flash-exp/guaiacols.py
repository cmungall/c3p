"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacol(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is a phenol with a methoxy group at the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise.
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for guaiacol substructure
    # Matches a 6-membered aromatic ring with an OH and a OCH3 adjacent to each other
    guaiacol_pattern = Chem.MolFromSmarts('c1([OH])c([OCH3])ccccc1')
    if guaiacol_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Check if the molecule has the guaiacol substructure
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Has Guaiacol substructure"
    else:
        return False, "Does not have Guaiacol substructure"
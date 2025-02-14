"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl group attached to a -CH2- group which is itself
    attached to an aromatic ring.

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

    # SMARTS pattern for an aromatic primary alcohol (aromatic-C-CH2-OH)
    aromatic_primary_alcohol_pattern = Chem.MolFromSmarts("c[CH2][OH]")

    # Find matches
    matches = mol.GetSubstructMatches(aromatic_primary_alcohol_pattern)

    if matches:
        return True, "Aromatic primary alcohol found."
    else:
       return False, "No aromatic primary alcohol found."
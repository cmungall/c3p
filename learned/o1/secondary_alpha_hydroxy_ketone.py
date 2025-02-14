"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a hydroxy group attached to the alpha carbon adjacent to a ketone group,
    where the alpha carbon is secondary (attached to two carbons and one hydrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for secondary alpha-hydroxy ketone
    pattern = Chem.MolFromSmarts('[CH1]([OH])[CH0](=O)')

    if mol.HasSubstructMatch(pattern):
        return True, "Contains secondary alpha-hydroxy ketone substructure"
    else:
        return False, "Does not contain secondary alpha-hydroxy ketone substructure"
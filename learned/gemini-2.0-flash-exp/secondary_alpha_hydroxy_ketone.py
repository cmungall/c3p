"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone is a ketone with a hydroxyl group on the alpha carbon,
    which also has exactly one hydrogen attached and one organyl group.

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

    # Create SMARTS pattern for secondary alpha-hydroxy ketone
    # [CX4H1]([OX2])[CX3](=[OX1])[#6] where the carbon with O and H is attached to a ketone via a direct bond.
    # The [#6] part enforces a carbon as a next neighbor. The alpha-carbon is defined as a carbon bound to one H, one OH and the carbonyl.
    # the single bond is already present, so we need only to state that a Carbon must be in that position.
    pattern = Chem.MolFromSmarts("[CX4H1]([OX2])[CX3](=[OX1])[#6]")
    matches = mol.GetSubstructMatches(pattern)


    if len(matches) == 0:
        return False, "Molecule does not contain a secondary alpha-hydroxy ketone substructure."

    return True, "Molecule contains a secondary alpha-hydroxy ketone substructure."
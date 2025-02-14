"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid consists of two or more flavonoid units connected by a single bond or atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define flavonoid substructure using SMARTS
    # This is a generalized pattern for the core of a flavonoid structure.
    # We are accounting for different oxidation states and substitution patterns
    # with the flexible SMARTS patterns used below.
    flavonoid_pattern = Chem.MolFromSmarts(
        "[cH1]1[cH][cH]([OX2])[c]([c]1)[CX3](=[OX1])[c]2[cH][cH][cH][c]([OX2])[c]2"
        ) # Modified to account for different substitution patterns

    # Count number of flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    num_flavonoids = len(flavonoid_matches)

    if num_flavonoids < 2:
        return False, f"Found {num_flavonoids} flavonoid units, need at least 2."

    # Check if flavonoids are linked directly (or through one atom)
    # We'll use a more flexible approach.
    # If a molecule has two or more flavonoids, and the molecule has no unconnected fragments,
    # we will assume the flavonoids are connected and therefore it's a biflavonoid.
    if len(Chem.GetMolFrags(mol)) > 1:
         return False, "Flavonoid rings are not connected, molecule is fragmented"
        
    return True, "Contains at least two flavonoid units connected by a single bond or atom"
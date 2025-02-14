"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to two carbon atoms and one hydrogen atom (implicitly).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary amine, where the nitrogen is bonded to 2 carbons
    # and the bond to the H is not explicit, meaning a single bond with one implicit hydrogen.
    # [NX3] is for a nitrogen with 3 bonds total, and the bond to hydrogen is not explicit.
    # [!#1] means non hydrogen, so not a nitrogen with a hydrogen
    # [CX4;!R] is a carbon with 4 single bonds, not in a ring.
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;!R]([CX4;!R])([CX4;!R])")

    # Check for at least one match
    if not mol.HasSubstructMatch(secondary_amine_pattern):
         return False, "No secondary amine found (N with two carbons)"

    # Check if the nitrogen is quaternized
    quaternary_nitrogen_pattern = Chem.MolFromSmarts("[NX4]")
    if mol.HasSubstructMatch(quaternary_nitrogen_pattern):
         return False, "Nitrogen is quaternized, not a secondary amine"

    return True, "Molecule contains a secondary amine group"
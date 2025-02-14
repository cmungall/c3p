"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: Flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone in which the hydrogen at position 3 of the heterocyclic ring
    is replaced by a hydroxy group (i.e., it's a 3-hydroxyflavone).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonol core SMARTS pattern with 3-hydroxy group
    flavonol_smarts = '[#6]1:[#6]:[#6](-c2ccc(O)cc2):[#8]c3cc(O)ccc13'  # Flavonol core with 3-OH
    flavonol_pattern = Chem.MolFromSmarts(flavonol_smarts)
    if flavonol_pattern is None:
        return False, "Invalid SMARTS pattern for flavonol core"

    # Check if molecule contains the flavonol core
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavonol core with hydroxyl at position 3 characteristic of flavonols"
    else:
        return False, "Does not contain the flavonol core structure with 3-hydroxyl group"
"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group in the beta position on the 3rd carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the steroid core with beta-hydroxyl at position 3.
    # [H] is a wildcard for implied H atoms
    # The [C@H] indicates that stereochemistry is specified at that center.
    # The pattern below has a general structure for the steroid core.
    # It specifies the beta configuration using [C@H] at position 3
    # The connected oxygen has a hydrogen implied.
    steroid_pattern = Chem.MolFromSmarts("[H][C@]12[C]([H])([C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C@H]1([H])[C]3([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])2[H])([H])O")

    # Check if the molecule matches the pattern
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Molecule does not match the steroid core with beta-hydroxyl at position 3"
        
    return True, "Molecule matches the criteria for a 3beta-hydroxy steroid"
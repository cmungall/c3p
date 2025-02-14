"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid has a steroid core with a beta-configured hydroxyl group at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid core as a SMARTS pattern, using X4 for C and specifying ring fusion
    steroid_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4]2[CX4]3[CX4]4[CX4]1[CX4]5[CX4]2[CX4]([CX4]3)[CX4]45")
    if not mol.HasSubstructMatch(steroid_core_pattern):
          return False, "No steroid core structure detected"

    # Define the 11-beta hydroxyl group SMARTS pattern, using the specific stereo bond
    # [C@H] indicates beta (up) config (it depends on the specific molecule, but in all examples, beta
    # is represented with this stereochemistry)
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[CX4]1[CX4]2[CX4]3[C@H](O)[CX4]4[CX4]1[CX4]5[CX4]2[CX4]([CX4]3)[CX4]45")
    
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11-beta-hydroxy group found"


    return True, "Molecule has a steroid core with an 11-beta-hydroxyl group"
"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for cyclohexenone
    # This pattern specifies:
    # - A six-membered ring with one double bond (C=C)
    # - A carbonyl group (C=O)
    # - Both are part of the same ring.
    cyclohexenone_pattern = Chem.MolFromSmarts("[C;R1](=[O])[C;R1]=[C;R1][C;R1][C;R1][C;R1]")
    
    
    # Check for match
    if not mol.HasSubstructMatch(cyclohexenone_pattern):
        return False, "Molecule does not contain a cyclohexenone ring"

    return True, "Molecule is a cyclohexenone"
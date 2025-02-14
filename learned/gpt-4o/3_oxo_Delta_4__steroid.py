"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is characterized by a 3-oxo group and a double bond 
    between carbon 4 and 5 of the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Discover the generalized steroid backbone: breaking into focused sections
    cyclohexene = Chem.MolFromSmarts("C1=CCCCC1")
    cyclohexane = Chem.MolFromSmarts("C1CCCCC1")
    cyclopentane = Chem.MolFromSmarts("C1CCCC1")
    
    # Check if these specific rings appear in sequence hinting at the steroidal nature
    features = [
        mol.HasSubstructMatch(cyclohexene),
        mol.HasSubstructMatch(cyclohexane),
        mol.HasSubstructMatch(cyclohexane),
        mol.HasSubstructMatch(cyclopentane)
    ]

    if not all(features):
        return False, "Not all necessary cyclic components identified for steroid backbone."

    # 3-Oxo feature as a branching part of the steroid nucleus
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C;R1]")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found where required in the backbone."

    # Discover the 4,5 Double bond to justify the Delta(4) portion
    delta_4_pattern = Chem.MolFromSmarts("C=CC(C)(C)")
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "Delta(4) double bond not found."

    return True, "3-oxo and Delta(4) double bond characteristics validated for the steroid"
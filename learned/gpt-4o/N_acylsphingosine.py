"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
"""

from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is composed of sphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingosine backbone pattern:
    # This is approximate, focusing on the chain with double bonds and hydroxyl and amine groups
    # The pattern includes a hydrocarbon chain, OH groups and NH
    sphingosine_pattern = Chem.MolFromSmarts("CO[C@@H](NC=O)CC\\C=C\\[C@H](O)CCCCC")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"

    # Look for amide bond pattern where nitrogen is connected to carbonyl carbon of a fatty acyl chain
    amide_pattern = Chem.MolFromSmarts("NC=O")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found connecting nitrogen to fatty acyl group"

    return True, "Contains sphingosine backbone with an amide bond connecting to a fatty acyl group"
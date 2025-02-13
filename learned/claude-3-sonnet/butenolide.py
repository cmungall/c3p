"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:30810 butenolide
A gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a butenolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-furanone skeleton
    furanone_pattern = Chem.MolFromSmarts("O=C1OCC=C1")
    if not mol.HasSubstructMatch(furanone_pattern):
        return False, "Missing 2-furanone skeleton"

    # Check for substitutions on the furanone ring
    substituted_pattern = Chem.MolFromSmarts("O=C1OC(*)C=C1*")
    if not mol.HasSubstructMatch(substituted_pattern):
        return True, "Unsubstituted 2-furanone skeleton (butenolide)"

    # Check for specific substitutions (optional)
    # ...

    return True, "Contains a substituted 2-furanone skeleton (butenolide)"
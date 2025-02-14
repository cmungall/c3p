"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:34831 thiol
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiol(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a thiol group (-SH) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thiol group (-SH) pattern
    thiol_pattern = Chem.MolFromSmarts("[SH]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)

    # Look for carbon-sulfur bond pattern
    c_s_pattern = Chem.MolFromSmarts("[C][S]")
    c_s_matches = mol.GetSubstructMatches(c_s_pattern)

    # A thiol must have at least one thiol group and one carbon-sulfur bond
    if len(thiol_matches) == 0 or len(c_s_matches) == 0:
        return False, "No thiol group or carbon-sulfur bond found"

    # Check if the thiol group is attached to a carbon atom
    for thiol_idx in thiol_matches:
        thiol_atom = mol.GetAtomWithIdx(thiol_idx)
        for neighbor in thiol_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                return True, "Contains a thiol group (-SH) attached to a carbon atom"

    return False, "Thiol group not attached to a carbon atom"
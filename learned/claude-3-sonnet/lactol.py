"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies chemical entities of the class lactol:
Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group. 
They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for lactol pattern: cyclic ether with OH group
    lactol_pattern = Chem.MolFromSmarts("[OH,OH2][C;r]1[O;r][C;r]=[C;r][C;r]=[C;r][C;r]1")
    matches = mol.GetSubstructMatches(lactol_pattern)
    if matches:
        return True, "Molecule contains a cyclic hemiacetal (lactol) structure"

    # Look for unsaturated lactol pattern
    unsaturated_lactol_pattern = Chem.MolFromSmarts("[OH,OH2][C;r]1[O;r][C;r]=[C;r][C;r]=[C;r]1")
    matches = mol.GetSubstructMatches(unsaturated_lactol_pattern)
    if matches:
        return True, "Molecule contains an unsaturated cyclic hemiacetal (lactol) structure"

    # If no matches found, it's not a lactol
    return False, "No lactol or unsaturated lactol structure found"
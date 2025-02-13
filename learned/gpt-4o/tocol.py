"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton substituted at position 2 
    by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chromanol core pattern
    chromanol_pattern = Chem.MolFromSmarts("Oc1cccc2c1CCCC2")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chromanol core found"
    
    # Check for substituent at position 2 (ignoring stereochemistry)
    # We assume indistinct pattern considering variable isoprenoid chains.
    substituent_pattern = Chem.MolFromSmarts("[C]c1cc(O)ccc1[C;!R](C)CCCC")
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "Position 2 is not properly substituted"

    # Check for a hydrocarbon chain pattern indicating isoprenoid units
    # Simplified pattern for isoprenoids: "(C)(C)C"
    isoprenoid_pattern = Chem.MolFromSmarts("[C](C)(C)C")
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_pattern)
    if len(isoprenoid_matches) < 3:
        return False, "Insufficient isoprenoid units"

    return True, "Contains chromanol core with appropriate substitution and isoprenoid units"
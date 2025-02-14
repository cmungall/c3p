"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of dimethylisoalloxazine with a substituent on the 10 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core isoalloxazine structure with variations on positions 7 and 8
    # The [NX3] is for an NH or a positive charged N at position 10
    core_pattern = Chem.MolFromSmarts("c1cc2[nX3]c3c([nX]c(=O)[nX]c3=O)n([!H])c2cc1")
    if not mol.HasSubstructMatch(core_pattern):
         return False, "Core dimethylisoalloxazine structure not found"

    # Define the isoalloxazine with 7 and 8 positions substituted with anything, and 
    # check for any substitution on the nitrogen at position 10.
    # this pattern does not enforce 7 and 8 to be methyl groups.
    
    full_pattern = Chem.MolFromSmarts("c1c([!H])c2[nX3]c3c([nX]c(=O)[nX]c3=O)n([!H])c2c([!H])c1")
    matches = mol.GetSubstructMatches(full_pattern)
    if len(matches) == 0:
        return False, "Core structure not found or not correctly substituted"

    # Get the nitrogen atom (N10)
    n10_match = mol.GetSubstructMatches(Chem.MolFromSmarts("[nX3]"))
    
    if not n10_match:
        return False, "N10 not found"

    for match in n10_match:
        n10_atom = mol.GetAtomWithIdx(match[0])
        if n10_atom.GetTotalDegree() > 2:
            return True, "Contains dimethylisoalloxazine core with a substituent on the 10 position"

    return False, "No substituent found at N10"
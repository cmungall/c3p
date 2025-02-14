"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group, -SH, is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alkanethiol.
    # This pattern looks for a sulfur atom bonded to a carbon atom. The carbon atom should have
    # 3 other connections, to C or H
    alkanethiol_pattern1 = Chem.MolFromSmarts("[SH][CX4]([#6,#1])([#6,#1])[#6,#1]") # covers the normal case of SH attached to CH2 or CH
    alkanethiol_pattern2 = Chem.MolFromSmarts("[SH][CH3]") # covers the case of a methyl group
    alkanethiol_pattern3 = Chem.MolFromSmarts("[SH][CX3]([#6,#1])=[#6]") # handles SH on a carbon adjacent to a C=C

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(alkanethiol_pattern1) or mol.HasSubstructMatch(alkanethiol_pattern2) or mol.HasSubstructMatch(alkanethiol_pattern3):
        
        matches = []
        if mol.HasSubstructMatch(alkanethiol_pattern1):
            matches.extend(mol.GetSubstructMatches(alkanethiol_pattern1))
        if mol.HasSubstructMatch(alkanethiol_pattern2):
            matches.extend(mol.GetSubstructMatches(alkanethiol_pattern2))
        if mol.HasSubstructMatch(alkanethiol_pattern3):
            matches.extend(mol.GetSubstructMatches(alkanethiol_pattern3))


        if matches:
           return True, "Contains a sulfanyl group (-SH) attached to an alkyl group."
            
    return False, "Does not contain a sulfanyl group (-SH) attached to an alkyl group."
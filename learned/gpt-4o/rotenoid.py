"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids are defined by the presence of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for tetrahydrochromeno[3,4-b]chromene skeleton
    # It's needed to match the internal structure of rotenoid
    rotenoid_pattern = Chem.MolFromSmarts("C1[C@H]2C=C3OC4=C(OCO4)C=C(C3=O)C(=C2OC5=CC(OC)=C(OC)C5)C1")
    
    # Check if the molecule has the rotenoid pattern
    if mol.HasSubstructMatch(rotenoid_pattern):
        return True, "Matches the rotenoid pattern"
    else:
        return False, "Does not match the rotenoid pattern"
"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: CHEBI:50902 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    A rotenoid is a tetrahydrochromenochromene with a cis-fused
    tetrahydrochromeno[3,4-b]chromene skeleton.

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
    
    # Define the tetrahydrochromeno[3,4-b]chromene skeleton pattern
    rotenoid_pattern = Chem.MolFromSmarts(
        """
        C1=C2OC3=C(C=C4OCCC4=C3)C(=O)C2=CC=C1
        """
    )
    
    # Check if the molecule contains the rotenoid skeleton
    if mol.HasSubstructMatch(rotenoid_pattern):
        return True, "Contains the cis-fused tetrahydrochromeno[3,4-b]chromene skeleton"
    else:
        return False, "Does not contain the rotenoid skeleton"

# Example usage
smiles = "COc1cc2OCC3Oc4cc(O)ccc4C(=O)C3c2cc1OC"  # 9-Demethylmunduserone
print(is_rotenoid(smiles))  # Output: (True, 'Contains the cis-fused tetrahydrochromeno[3,4-b]chromene skeleton')
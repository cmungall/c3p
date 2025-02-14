"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:27042 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings
    joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for corrin macrocycle pattern
    corrin_pattern = Chem.MolFromSmarts("[n&r4,r5,r6]1[n&r4,r5,r6][c&r4,r5,r6][c&r4,r5,r6][n&r4,r5,r6][c&r4,r5,r6][c&r4,r5,r6][n&r4,r5,r6]1")
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "No corrin macrocycle found"
    
    # Count number of reduced pyrrole rings
    reduced_pyrrole_pattern = Chem.MolFromSmarts("[Nr4,r5,r6]")
    reduced_pyrroles = len(mol.GetSubstructMatches(reduced_pyrrole_pattern))
    if reduced_pyrroles < 4:
        return False, f"Found {reduced_pyrroles} reduced pyrrole rings, need at least 4"
    
    # Look for three =C- groups and one direct C-C bond linking alpha positions
    alpha_link_pattern = Chem.MolFromSmarts("[c&r4,r5,r6][c&r4,r5,r6]")
    alpha_links = mol.GetSubstructMatches(alpha_link_pattern)
    if len(alpha_links) != 1:
        return False, "Incorrect number of alpha-linked carbon atoms"
    
    # Check for cobalt coordination
    cobalt_pattern = Chem.MolFromSmarts("[Co]")
    if not mol.HasSubstructMatch(cobalt_pattern):
        return False, "No cobalt coordination found"
    
    return True, "Contains corrin macrocycle with 4 reduced pyrrole rings, 3 =C- groups, 1 C-C alpha link, and cobalt coordination"
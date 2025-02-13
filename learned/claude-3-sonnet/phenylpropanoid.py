"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:37134 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring
    if not mol.GetAromaticRings():
        return False, "No aromatic rings found"

    # Check for phenylpropane skeleton
    phenylpropane_pattern = Chem.MolFromSmarts("[c:1]1[c;r6]([c;r6][c;r6][c:2]1)[C:3]=[C:4][C:5]")
    match = mol.GetSubstructMatches(phenylpropane_pattern)
    if not match:
        return False, "No phenylpropane skeleton found"

    # Check for additional rings, oxy substitutions, etc.
    ring_patterns = [Chem.MolFromSmarts("[Or:6]"), Chem.MolFromSmarts("[c:6]1[c;r6][c;r6][c;r6][c;r6][c;r6]1")]
    extra_features = any(mol.HasSubstructMatch(pattern) for pattern in ring_patterns)

    # Basic phenylpropanoids have a phenylpropane skeleton and aromatic ring
    if not extra_features:
        return True, "Contains phenylpropane skeleton and aromatic ring"

    # More complex phenylpropanoids may have additional rings, oxy substitutions, etc.
    return True, "Contains phenylpropane skeleton, aromatic ring, and additional features"
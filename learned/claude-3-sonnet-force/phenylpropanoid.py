"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26266 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

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
    
    # Look for phenylpropane backbone (benzene ring with 3-carbon chain)
    backbone_pattern = Chem.MolFromSmarts("[c]1ccccc1CCC")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No phenylpropane backbone found"
    
    # Check for aromatic ring
    if sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) < 6:
        return False, "Not sufficiently aromatic"
    
    # Check for common substituents (hydroxy, methoxy, ester, etc.)
    substituents = ['O', 'OC', 'C(=O)O', 'C(=O)OC', 'c1ccccc1']
    substituent_pattern = Chem.MolFromSmarts("|".join(["[$({})]/[#6]".format(sub) for sub in substituents]))
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "No common phenylpropanoid substituents found"
    
    return True, "Contains phenylpropane backbone with aromatic ring and common substituents"
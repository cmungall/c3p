"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:37335 xanthophyll
A subclass of carotenoids consisting of the oxygenated carotenes.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid containing oxygen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carotenoid backbone (long chain of conjugated double bonds)
    carotenoid_pattern = Chem.MolFromSmarts("[C;R]=[C;R][C;R]=[C;R][C;R]=[C;R][C;R]=[C;R][C;R]=[C;R]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
        return False, "No carotenoid backbone found"
    
    # Check for oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chain too short for a carotenoid"

    # Check molecular weight - xanthophylls typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for xanthophyll"

    return True, "Contains a carotenoid backbone with oxygen atoms"
"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more hydroxyl groups replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for 5 or 6 membered rings with C and O atoms
    sugar_ring_pattern = Chem.MolFromSmarts("[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]1")
    sugar_ring_pattern2 = Chem.MolFromSmarts("[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3][CX4,CX3]1")
    if not mol.HasSubstructMatch(sugar_ring_pattern) and not mol.HasSubstructMatch(sugar_ring_pattern2):
        return False, "No sugar ring structure found"
        
    # Look for any N atom bonded to a ring carbon
    amino_group_pattern = Chem.MolFromSmarts("[NX3;!H0]-[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]1")
    amino_group_pattern2 = Chem.MolFromSmarts("[NX3;!H0]-[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3][CX4,CX3]1")
    
    amino_matches = mol.GetSubstructMatches(amino_group_pattern)
    amino_matches2 = mol.GetSubstructMatches(amino_group_pattern2)
    if len(amino_matches) == 0 and len(amino_matches2) == 0:
       return False, "No amino group substitution found"

    # If both conditions are met, it is an amino sugar
    return True, "Contains a sugar ring with one or more amino group substitutions."
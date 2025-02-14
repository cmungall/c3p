"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:35349 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as any compound based on the cyclopenta[a]phenanthrene carbon skeleton,
    with methyl groups at C-10 and C-13, and often an alkyl group at C-17. Ring expansions,
    contractions, and bond scissions are allowed.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopenta[a]phenanthrene backbone
    backbone_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)CC2)C)CC1"
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No cyclopenta[a]phenanthrene backbone found"
    
    # Look for methyl groups at C-10 and C-13
    c10_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)CC2)C(C)CC1")
    c13_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)CC2)CC(C)C1")
    if not mol.HasSubstructMatch(c10_pattern) or not mol.HasSubstructMatch(c13_pattern):
        return False, "Missing methyl groups at C-10 and/or C-13"
    
    # Check for alkyl group at C-17
    c17_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C[C@@H]([C@H](C)CC)C)CC2)CC1")
    c17_match = mol.HasSubstructMatch(c17_pattern)
    
    # Allow for ring expansions, contractions, and bond scissions
    ring_expansion_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)[C@@H]2CC1)[C@@H]5[C@@H]5")
    ring_contraction_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)[C@@H]2CC1)[C@@H]5C5")
    bond_scission_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)[C@@H]2C[C@@H]1[C@@H]5[C@@H]5")
    ring_expansion_match = mol.HasSubstructMatch(ring_expansion_pattern)
    ring_contraction_match = mol.HasSubstructMatch(ring_contraction_pattern)
    bond_scission_match = mol.HasSubstructMatch(bond_scission_pattern)
    
    if c17_match and (ring_expansion_match or ring_contraction_match or bond_scission_match):
        return True, "Contains cyclopenta[a]phenanthrene backbone with methyl groups at C-10 and C-13, and alkyl group at C-17, with possible ring expansions, contractions, or bond scissions"
    elif c17_match:
        return True, "Contains cyclopenta[a]phenanthrene backbone with methyl groups at C-10 and C-13, and alkyl group at C-17"
    else:
        return True, "Contains cyclopenta[a]phenanthrene backbone with methyl groups at C-10 and C-13"
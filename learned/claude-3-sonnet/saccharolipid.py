"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI:36975 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid that contains a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Carbohydrate patterns
    monosaccharide_patterns = ["OC[C@H](O)[C@H](O)[C@H](O)CO", # aldose
                               "OC[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", # furanose
                               "OC[C@@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@H]1O"] # pyranose
    
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2]C[OX2]")
    
    # Look for carbohydrate moieties
    has_carbohydrate = False
    for pattern in monosaccharide_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_carbohydrate = True
            break
    
    if not has_carbohydrate:
        return False, "No carbohydrate moiety found"
    
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bonds found"
    
    # Lipid patterns
    linear_lipid_pattern = Chem.MolFromSmarts("[C;H3][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2]")
    branched_lipid_pattern = Chem.MolFromSmarts("[C;H3][C;H2][C;H2]([C;H3])[C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2]")
    cyclic_lipid_pattern = Chem.MolFromSmarts("[C;H2]1[C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2]1")
    
    # Look for lipid moieties
    has_lipid = mol.HasSubstructMatch(linear_lipid_pattern) or \
                mol.HasSubstructMatch(branched_lipid_pattern) or \
                mol.HasSubstructMatch(cyclic_lipid_pattern)
    
    if not has_lipid:
        return False, "No lipid moiety found"
    
    return True, "Contains both carbohydrate and lipid moieties"
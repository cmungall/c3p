"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:27327 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids with a polyene backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check number of carbons: should be at least 30
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 30:
        return False, "Too few carbon atoms for a xanthophyll."

    # Check for a longer polyene backbone, common for carotenoids
    # More permissive pattern for conjugated double bonds
    carotenoid_pattern = Chem.MolFromSmarts("[CX3]=[CX3]~[CX3]=[CX3]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
      return False, "No basic carotenoid-like backbone structure found."

    # Check for oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms present, not a xanthophyll."

    # Check for oxygen attached to sp3 carbon or double bonded to C
    oxygenated_carbon_pattern = Chem.MolFromSmarts("[CX4]-[OX2]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    oxygenated_groups = mol.GetSubstructMatches(oxygenated_carbon_pattern) + mol.GetSubstructMatches(carbonyl_pattern)


    if len(oxygenated_groups) < 1:
         return False, "No characteristic oxygen-containing functional groups found."


    # Check for glycoside or sulfates using a more permissive SMARTS
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H]([OX2])[C@@H](CO)[OX1]")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)([OX1])([OX1])")
    if mol.HasSubstructMatch(glycoside_pattern) or mol.HasSubstructMatch(sulfate_pattern):
      return True, "Contains a carotenoid-like backbone with characteristic oxygen-containing functional groups including glycosides or sulfates."
    
    
    # Check molecular weight - xanthophylls are typically > 400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a xanthophyll"

    return True, "Contains a carotenoid-like backbone with oxygen-containing functional groups."
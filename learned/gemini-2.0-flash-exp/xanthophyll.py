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
    
    # Basic Carotenoid backbone - conjugated double bonds with methyl groups
    carotenoid_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
        return False, "No carotenoid-like backbone structure found."

    # Check for the presence of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
      return False, "No oxygen atoms, so cannot be a xanthophyll."

    # check number of carbons, should be roughly 40
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 30 or carbon_count > 50:
        return False, "Number of carbon atoms not within typical range for xanthophylls."

    #Check for functional groups which contain oxygen: hydroxyl, carbonyl, epoxides, and glycosides
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    carbonyl_pattern = Chem.MolFromSmarts("C=[OX1]")
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    glycoside_pattern = Chem.MolFromSmarts("OC[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1") # Simplified glycoside pattern

    if not (mol.HasSubstructMatch(hydroxy_pattern) or mol.HasSubstructMatch(carbonyl_pattern) or mol.HasSubstructMatch(epoxide_pattern) or mol.HasSubstructMatch(glycoside_pattern)):
      return False, "No oxygen-containing functional groups found which are typical for xanthophylls."

    return True, "Contains a carotenoid-like backbone with oxygen-containing functional groups."
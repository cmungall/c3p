"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid that contains a carbohydrate moiety.

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

    # Define a pattern for carbohydrate moiety (sugar rings)
    # Sugar rings are typically five or six-membered rings with oxygen and multiple hydroxyl groups
    sugar_pattern = Chem.MolFromSmarts("[OX2H][CX4H1,O][CX4H1,O][CX4H1,O][CX4H1,O][OX2H]")  # Six-membered ring
    five_ring_sugar_pattern = Chem.MolFromSmarts("[OX2H][CX4H1,O][CX4H1,O][CX4H1,O][OX2H]")  # Five-membered ring
    
    # Check for presence of sugar moiety
    has_sugar = mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(five_ring_sugar_pattern)
    if not has_sugar:
        return False, "No carbohydrate (sugar) moiety found"

    # Define a pattern for lipid moiety (long hydrocarbon chains)
    # A simple pattern is a chain of at least 8 carbons (aliphatic chain)
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCC")  # Eight or more consecutive carbons

    # Check for presence of lipid moiety
    has_lipid = mol.HasSubstructMatch(lipid_pattern)
    if not has_lipid:
        return False, "No lipid moiety found (long hydrocarbon chains)"

    # Additional checks could involve verifying ester or amide bonds linking the lipid and sugar
    ester_or_amide_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,C,N]")  # Ester or amide bond
    has_linkage = mol.HasSubstructMatch(ester_or_amide_pattern)
    if not has_linkage:
        return False, "No ester or amide linkage found between lipid and sugar moieties"

    return True, "Contains both carbohydrate and lipid moieties linked via ester or amide bonds"
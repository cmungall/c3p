"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid contains a carbohydrate moiety linked to a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns to search for: carbohydrate (sugar) moiety and lipid (fatty acid) chains
    
    # A simple carbohydrate pattern: a cyclic structure with multiple oxygens
    carbohydrate_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1")  # generic hexose-like structure
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No recognizable carbohydrate moiety found"

    # A simple fatty acid or lipid chain: long hydrocarbon chain possibly terminating in an ester or amide linkage
    lipid_chain_pattern = Chem.MolFromSmarts("C(=O)[O,N][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")  # Fatty acid ester/amide
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No recognizable lipid/fatty acid component found"

    # Check if both carbohydrate and lipid are part of the same molecule
    carb_matches = mol.GetSubstructMatches(carbohydrate_pattern)
    lipid_matches = mol.GetSubstructMatches(lipid_chain_pattern)
    
    # A simplistic check assuming that both matches exist and bonds are present between them
    for carb_match in carb_matches:
        for lipid_match in lipid_matches:
            if any(bond.IsInRing() for bond in mol.GetBondsBetweenAtoms(carb_match[0], lipid_match[0])):
                return True, "Contains both carbohydrate and lipid components correctly linked"
    
    return False, "Carbohydrate and lipid components not sufficiently linked in structure"

# This is a rough implementation for model classification. The patterns and linkages need refinement for robust applications.
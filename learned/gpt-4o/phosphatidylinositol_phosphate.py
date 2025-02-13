"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for inositol ring pattern (six-membered ring with hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C(O)1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring structure found"

    # Check for one or more phosphate groups (-P(=O)(O)O-)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # Confirm presence of long fatty acid chains (at least two of significant length)
    long_chain_pattern = Chem.MolFromSmarts("C(=O)OC(C)C")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Not enough long carbon chains"

    # Check molecular weight which is typically large for phosphatidylinositol phosphates
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:  # arbitrary threshold for large molecules in this class
        return False, f"Molecular weight too low for phosphatidylinositol phosphate: {mol_wt}"

    return True, "Contains glycerol backbone, inositol ring, phosphate groups, and long carbon chains typical of phosphatidylinositol phosphates"
"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule could potentially be classified as a semisynthetic derivative
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Consider known modification groups and combinations commonly used in semisynthetic derivatives
    # More focused structures than the previous attempt.
    
    # Check for specific complex ester patterns, unusual in natural products
    complex_ester_group = Chem.MolFromSmarts("COC(=O)C")
    if mol.HasSubstructMatch(complex_ester_group):
        return True, "Contains complex ester group, indicative of chemical modification"

    # Combined ether and ester groups close in structure
    ether_ester_combination = Chem.MolFromSmarts("C(OC)C(=O)O")
    if mol.HasSubstructMatch(ether_ester_combination):
        return True, "Contains ether and ester combination, indicative of semisynthetic processing"

    # Include sugar-like structures attached to non-sugar core (common indication of derivatization)
    sugar_derivative_pattern = Chem.MolFromSmarts("[C@H]1([O][C@H]([C@@H](O)[C@H]1O)O[C@H]1C[C@@H](O)[C@@H]1O)C")
    if mol.HasSubstructMatch(sugar_derivative_pattern):
        return True, "Contains sugar derivative, indicative of semisynthetic derivation"

    # Include larger macrocyclic esters or lactone rings that are rare in natural products but common in derivatization
    macrocyclic_lactone = Chem.MolFromSmarts("C1(OC(=O)C)OC1")
    if mol.HasSubstructMatch(macrocyclic_lactone):
        return True, "Contains macrocyclic lactone, indicative of synthetic derivation"

    # Mixed amide and ester groups
    amide_ester_combination = Chem.MolFromSmarts("C(=O)NC(=O)O")
    if mol.HasSubstructMatch(amide_ester_combination):
        return True, "Contains amide and ester combination, potential semisynthetic derivative"

    # Prioritize known complex synthetic modifications such as acetals in complex structures
    acetal_like_structure = Chem.MolFromSmarts("C(OC)(OC)C")
    if mol.HasSubstructMatch(acetal_like_structure):
        return True, "Contains acetal-like structure, indicates synthetic methodology"

    return False, "Molecule does not have typical refined indicators of semisynthetic derivation"
"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
"""
from rdkit import Chem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of a central C20 backbone
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Contains {carbon_count} carbons; requires at least 20 carbons for backbone"
    
    # Check for presence of key functional groups: epoxide, multiple hydroxyls, and carboxylic acid
    # Epoxide pattern with allowance for asymmetric structures if needed
    epoxide_pattern = Chem.MolFromSmarts("[C;R1]1[O;R1][C;R1]1")  
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")  # Carboxylic acid group
    
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)

    # Ensure it has at least one epoxide, at least two hydroxyls, and one carboxylic acid
    if not has_epoxide or hydroxyl_count < 2 or not has_carboxylic_acid:
        return False, "Must contain epoxide group, multiple hydroxyls, and a carboxylic acid group"
    
    # Avoid classic icosanoid motifs
    prostanoid_pattern = Chem.MolFromSmarts("C1=CCCCC1C(=O)O")  # Generalized prostanoids ring
    leukotriene_pattern = Chem.MolFromSmarts("CC(C)(C)C[CX4]")  # Branching pattern specific

    if mol.HasSubstructMatch(prostanoid_pattern) or mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Substructure resembles classic icosanoids (prostanoid or leukotriene patterns detected)"
    
    # If all criteria are met, classify as nonclassic icosanoid
    return True, "Classified as nonclassic icosanoid with adequate C20 backbone and functionalization"

# Note: This classification is based on identifying structural elements typically 
# found in nonclassic icosanoids and may not capture all derivatives.
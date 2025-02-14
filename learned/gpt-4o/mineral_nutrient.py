"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic compound essential for the human body,
    typically consisting of an inorganic anion and a metal cation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if classified as a mineral nutrient;
                   False and reason otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define valid SMARTS patterns for metal cations
    metal_cations_patterns = [
        Chem.MolFromSmarts("[K+]"),
        Chem.MolFromSmarts("[Na+]"),
        Chem.MolFromSmarts("[Ca+2]"),
        Chem.MolFromSmarts("[Mg+2]"),
        Chem.MolFromSmarts("[Fe+3]"),
        Chem.MolFromSmarts("[Zn+2]"),
        Chem.MolFromSmarts("[Al+3]"),
        Chem.MolFromSmarts("[Ba+2]"),
        Chem.MolFromSmarts("[Cs+]"),
        Chem.MolFromSmarts("[La+3]"),
        Chem.MolFromSmarts("[Sb+5]")
    ]
    
    # Define valid SMARTS patterns for common anionic groups
    anions_patterns = [
        Chem.MolFromSmarts("[Cl-]"),
        Chem.MolFromSmarts("[F-]"),
        Chem.MolFromSmarts("[O-]S(=O)(=O)[O-]"),  # sulfate
        Chem.MolFromSmarts("[O-]P(=O)([O-])[O-]"),  # phosphate
        Chem.MolFromSmarts("[O-]C(=O)[O-]"),  # carbonate
        Chem.MolFromSmarts("[O-]N(=O)[O-]"),  # nitrate
        Chem.MolFromSmarts("[Si]([O-])([O-])([O-])[O-]")  # silicate
    ]
    
    # Check for presence of at least one metal cation
    cation_found = any(mol.HasSubstructMatch(pattern) for pattern in metal_cations_patterns)
    if not cation_found:
        return False, "No essential metal cations found"
    
    # Check for presence of at least one anion or inorganic group
    anion_found = any(mol.HasSubstructMatch(pattern) for pattern in anions_patterns)
    if not anion_found:
        return False, "No recognizable anionic groups found"
    
    # If both components are present, classify as a potential mineral nutrient
    return True, "Contains essential metal cation and anionic counterpart typical of mineral nutrients"
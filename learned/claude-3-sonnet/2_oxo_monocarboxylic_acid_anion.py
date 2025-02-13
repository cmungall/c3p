"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:38478 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    Must have a carboxylate group with an oxo group at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate group pattern ([O-]C(=O)-)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if not carboxylate_matches:
        return False, "No carboxylate group found"
    
    if len(carboxylate_matches) > 1:
        return False, "More than one carboxylate group found"

    # Look for the 2-oxo pattern: [O-]C(=O)C(=O)R
    oxo_acid_pattern = Chem.MolFromSmarts("[O-]C(=O)C(=O)[#6,#7,#8,#9,#16,#15]")
    if not mol.HasSubstructMatch(oxo_acid_pattern):
        return False, "No 2-oxo group adjacent to carboxylate found"
    
    # Make sure there's only one ketone group
    ketone_pattern = Chem.MolFromSmarts("[#6]C(=O)[#6]")
    additional_ketones = mol.GetSubstructMatches(ketone_pattern)
    
    # The 2-oxo group will be counted in additional_ketones, so we check if there are more than 1
    if len(additional_ketones) > 1:
        return False, "Additional ketone groups found beyond the 2-position"
    
    # Additional check to ensure the oxo group is at position 2
    # First get the carboxylate carbon atom
    carboxylate_carbon = carboxylate_matches[0][1]  # Index 1 is the carbon of C(=O)[O-]
    
    # Get all ketone carbons
    ketone_carbons = [match[1] for match in mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)"))]
    
    # Check if any ketone carbon is directly bonded to the carboxylate carbon
    for ketone_carbon in ketone_carbons:
        bond = mol.GetBondBetweenAtoms(carboxylate_carbon, ketone_carbon)
        if bond is not None:
            return True, "Contains carboxylate with oxo group at 2-position"
    
    return False, "Oxo group not at 2-position"
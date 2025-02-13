"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are characterized by a 4-membered beta-lactam ring fused to a 5-membered ring 
    containing sulfur, with substitutions at positions 3, 4, and 6.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a beta-lactam ring, a 4-membered ring containing nitrogen and carbonyl groups
    beta_lactam_pattern = Chem.MolFromSmarts('[NX3]1C(=O)C[CX3]1')
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check for the presence of a 5-membered ring containing sulfur and fused to the beta-lactam
    carbapenem_pattern = Chem.MolFromSmarts('C1OC(=O)C2SC=C(N2)C=C1')
    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "No adjacent sulfur-containing 5-membered ring fused to beta-lactam found"
    
    # Optionally, we can search for common substitutions around the rings if further precision is needed

    return True, "Contains the characteristic structure of a carbapenem antibiotic"
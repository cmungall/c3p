"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C][C][C]3[C]2[C][C][C]4[C][C][C][C][C]34")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # The pattern looks for a carbonyl (C=O) group in the A ring of the steroid
    oxo_pattern = Chem.MolFromSmarts("[CH2][CH2][C](=[O])[CH2][C](@[H])")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-alpha configuration
    # In 5-alpha steroids, the hydrogen at position 5 is below the plane (alpha)
    # This means the connection between rings A/B is trans
    # Looking for the specific connectivity pattern around C5
    alpha_5_pattern = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@@H]([C]2)[CH2][CH2]")
    if not mol.HasSubstructMatch(alpha_5_pattern):
        return False, "No 5-alpha configuration found"

    # Additional check for reasonable molecular weight range of steroids
    mol_wt = sum([atom.GetMass() for atom in mol.GetAtoms()])
    if mol_wt < 250 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for steroids"

    # Count carbons (steroids typically have 19+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19:
        return False, "Too few carbons for a steroid structure"

    return True, "Molecule contains 3-oxo-5alpha-steroid structure"
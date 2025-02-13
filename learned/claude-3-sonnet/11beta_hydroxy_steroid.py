"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11beta-hydroxy steroids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (4 fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 11-beta hydroxy group
    # The SMARTS pattern looks for:
    # - The specific carbon at position 11 in the steroid skeleton
    # - An OH group in beta configuration (specified by '@H')
    # Note: The exact stereochemistry is critical here
    hydroxy_11beta = Chem.MolFromSmarts("[C]~1~[C]~[C]~[C]~2~[C]~[C]~[C]~[C]~3~[C]~[C@@H](O)~[C]~[C]~4~[C]~[C]~[C]~[C]~4~[C]~3~[C]~2~1")
    
    if not mol.HasSubstructMatch(hydroxy_11beta):
        return False, "No 11-beta hydroxy group found"
    
    # Additional check for common substituents often found in 11beta-hydroxy steroids
    # Look for common functional groups like ketones, other hydroxyls, etc.
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    other_hydroxy = Chem.MolFromSmarts("[CH]O")
    
    if not (mol.HasSubstructMatch(ketone_pattern) or mol.HasSubstructMatch(other_hydroxy)):
        return False, "Missing typical steroid substituents"
        
    # Count carbons to ensure reasonable size for a steroid
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19 or carbon_count > 30:  # Most steroids have 19-30 carbons
        return False, f"Carbon count ({carbon_count}) outside typical steroid range"
    
    # Check for reasonable number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Too few rings for steroid structure"
        
    return True, "Contains steroid core with 11-beta hydroxy group and appropriate substituents"
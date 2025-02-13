"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
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

    # Check for 17-OH with alpha stereochemistry
    # [C] is carbon 17, [OH1] is the hydroxyl group, '@' indicates stereochemistry
    # The exact SMARTS pattern depends on the numbering convention used in the structure
    oh_17_alpha_pattern = Chem.MolFromSmarts('[C;R1]-[C;R1]([OH1])')
    
    if not mol.HasSubstructMatch(oh_17_alpha_pattern):
        return False, "No hydroxyl group at C17 position found"
    
    # Additional checks to verify it's a steroid structure:
    # Count rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"
    
    # Basic size check - steroids typically have at least 19 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 19:
        return False, "Too few carbons for steroid structure"
    
    # Check for sp3 hybridized carbons typical in steroid core
    sp3_carbons = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3)
    if sp3_carbons < 10:
        return False, "Insufficient sp3 carbons for steroid structure"

    return True, "Contains steroid core with 17-alpha hydroxyl group"
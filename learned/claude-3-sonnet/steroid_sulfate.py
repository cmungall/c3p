"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for steroid core (four fused rings - three 6-membered and one 5-membered)
    # This SMARTS pattern represents the basic steroid skeleton
    steroid_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core structure found"
    
    # Look for sulfate group (-OS(=O)(=O)O or -OS(=O)(=O)[O-])
    sulfate_pattern1 = Chem.MolFromSmarts("OS(=O)(=O)O")
    sulfate_pattern2 = Chem.MolFromSmarts("OS(=O)(=O)[O-]")
    
    has_sulfate1 = mol.HasSubstructMatch(sulfate_pattern1)
    has_sulfate2 = mol.HasSubstructMatch(sulfate_pattern2)
    
    if not (has_sulfate1 or has_sulfate2):
        return False, "No sulfate group found"
    
    # Count sulfate groups
    sulfate_matches1 = len(mol.GetSubstructMatches(sulfate_pattern1))
    sulfate_matches2 = len(mol.GetSubstructMatches(sulfate_pattern2))
    total_sulfates = sulfate_matches1 + sulfate_matches2
    
    # Verify that sulfate is attached to steroid core
    # Look for O-S bond where O is attached to carbon (part of steroid)
    sulfate_linkage = Chem.MolFromSmarts("[#6]-[#8]-[#16](=[#8])(=[#8])-[#8]")
    if not mol.HasSubstructMatch(sulfate_linkage):
        return False, "Sulfate group not properly connected to steroid core"
    
    # Count carbons to ensure reasonable size for steroid
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19:  # Most steroids have 19+ carbons
        return False, "Too few carbons for steroid structure"
        
    if total_sulfates == 1:
        return True, "Found steroid core with one sulfate group"
    else:
        return True, f"Found steroid core with {total_sulfates} sulfate groups"
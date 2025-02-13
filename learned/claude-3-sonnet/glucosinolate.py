"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: glucosinolate compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the core glucosinolate structure:
    # - Sulfur linked to sugar and C=N
    # - C=N double bond
    # - N-O-SO3 group
    core_pattern = Chem.MolFromSmarts('[#6]-[#16]-[#6]=[#7]-[#8]-[#16](=[#8])(=[#8])[#8]')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core glucosinolate structure (S-C=N-O-SO3)"

    # Check for pyranose sugar (glucose) pattern
    sugar_pattern = Chem.MolFromSmarts('[#6]1-[#8]-[#6]-[#6]-[#6]-[#6]1')
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "Missing pyranose sugar ring"
    
    # Count hydroxyl groups on sugar (should be multiple OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts('[#8H1]-[#6]')
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 4:  # Glucose has multiple OH groups
        return False, "Insufficient hydroxyl groups for glucose moiety"

    # Verify S-glycosidic linkage
    glycosidic_pattern = Chem.MolFromSmarts('[#6]1-[#8]-[#6]-[#6]-[#6]-[#6]1-[#16]')
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Missing S-glycosidic linkage"
    
    # Check for sulfate group
    sulfate_pattern = Chem.MolFromSmarts('[#16](=[#8])(=[#8])-[#8]')
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "Missing sulfate group"

    # Count key atoms to ensure reasonable composition
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if s_count < 2:  # Need at least 2 sulfur atoms
        return False, "Insufficient sulfur atoms"
    if o_count < 7:  # Need multiple oxygens for sugar + sulfate
        return False, "Insufficient oxygen atoms"
    if n_count < 1:  # Need 1 nitrogen
        return False, "Missing nitrogen atom"

    return True, "Contains glucosinolate core structure with glucose moiety and characteristic S-C=N-O-SO3 group"
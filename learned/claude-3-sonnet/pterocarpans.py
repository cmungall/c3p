"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: pterocarpans
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a characteristic 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic pterocarpan skeleton SMARTS pattern
    # Represents the core 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene structure
    pterocarpan_core = Chem.MolFromSmarts(
        "[#6]1-[#6]2-[#6]-[#6]3=,:[#6]-[#6](=,:[#6]-[#6]=,:[#6]3-[#8]-[#6]2)-[#8]-[#6]1"
    )
    
    if not mol.HasSubstructMatch(pterocarpan_core):
        return False, "Missing characteristic pterocarpan core structure"

    # Check for the presence of two benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = len(mol.GetSubstructMatches(benzene_pattern))
    if benzene_matches < 2:
        return False, "Missing required benzene rings"

    # Check for the characteristic oxygen bridges (furan and pyran rings)
    oxygen_bridge_pattern = Chem.MolFromSmarts("[#6]~[#8]~[#6]")
    oxygen_bridges = len(mol.GetSubstructMatches(oxygen_bridge_pattern))
    if oxygen_bridges < 2:
        return False, "Missing characteristic oxygen bridges"

    # Check for sp3 carbons at 6a and 11a positions
    sp3_carbon_pattern = Chem.MolFromSmarts("[CX4]")
    sp3_carbons = len(mol.GetSubstructMatches(sp3_carbon_pattern))
    if sp3_carbons < 2:
        return False, "Missing required sp3 carbons at 6a and 11a positions"

    # Verify the molecule is not a coumestan (which would have a C=O group)
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl group (appears to be a coumestan rather than pterocarpan)"

    # Additional check for the complete ring system
    tetracyclic_system = Chem.MolFromSmarts(
        "[#6]1-2-[#6]-[#6]3=,:[#6]-[#6](=,:[#6]-[#6]=,:[#6]3-[#8]-[#6]-2)-[#8]-[#6]-1"
    )
    if not mol.HasSubstructMatch(tetracyclic_system):
        return False, "Missing complete tetracyclic ring system"

    return True, "Contains characteristic 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton"
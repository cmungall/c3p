"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with specific structural characteristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count Rings (must be 4 or more)
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, "Less than 4 rings"
    
    # Count Carbons (must be about 30)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27 or c_count > 33:
         return False, f"Carbon count {c_count} is out of typical range (27-33)"
    
    # Count Oxygens (must be at least 4)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
         return False, f"Oxygen count {o_count} is too low"

    # Check for tetracyclic core (this is a simplified SMARTS, as it's difficult to make a general one)
    # This pattern just looks for a 4 ring fused system, without any constraints on the atoms in it
    core_pattern1 = Chem.MolFromSmarts("C12C3C4C1C5C2C3C45")
    core_pattern2 = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]1")
    
    if not mol.HasSubstructMatch(core_pattern1) and not mol.HasSubstructMatch(core_pattern2):
        return False, "Tetracyclic core not found"

    # Check for at least 4 oxygens on the ring system (more than just alcohol groups)
    # Check for common side chain patterns
    # This includes alcohols and carbonyls attached to the core rings
    oxygen_pattern = Chem.MolFromSmarts("[CX4]([OX2])[CX4]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    if (len(oxygen_matches) + len(carbonyl_matches)) < 3:
       return False, "Not enough oxygens on the core rings"

    # Check for possible isopropyl group at C20
    isopropyl_pattern = Chem.MolFromSmarts("[CX4](C)(C)")
    isopropyl_matches = mol.GetSubstructMatches(isopropyl_pattern)
    if len(isopropyl_matches) < 1 : # many cucurbitacins have the isopropyl, some don't.
          pass # don't fail, as isopropyl isn't mandatory, and may be substituted
    

    # If all checks passed, we consider it a cucurbitacin
    return True, "Likely a cucurbitacin based on structure features"
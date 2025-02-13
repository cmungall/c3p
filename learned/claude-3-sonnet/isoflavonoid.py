"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:24474 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core benzopyran structure with aryl at position 3
    # [#6] represents any carbon atom
    # The numbers in the SMARTS pattern indicate atom mapping for reference
    isoflavonoid_core = Chem.MolFromSmarts(
        """
        [#6]1:,[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]2-[#8]-[#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1-2
        """
    )
    
    # Check for core benzopyran structure
    if not mol.HasSubstructMatch(isoflavonoid_core):
        return False, "No benzopyran core structure found"

    # Look for aromatic ring attached at position 3
    # This pattern specifically looks for the isoflavonoid arrangement
    isoflavonoid_pattern = Chem.MolFromSmarts(
        """
        [#6]1:,[#6](=[O,CH2])-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#8]-[#6]:1-[#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1
        """
    )
    
    if not mol.HasSubstructMatch(isoflavonoid_pattern):
        return False, "No aryl substituent at position 3 of benzopyran"

    # Additional check to avoid false positives
    # Count number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for isoflavonoid structure"

    # Most isoflavonoids have oxygen-containing substituents
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Insufficient oxygen atoms for typical isoflavonoid"

    return True, "Contains benzopyran core with aryl substituent at position 3"
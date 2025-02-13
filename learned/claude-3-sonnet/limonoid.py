"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:64655 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a prototypical 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for furan ring
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"
    
    # Look for triterpene skeleton (4 rings and a high degree of oxygenation)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Less than 4 rings found"
    
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Less than 5 oxygens found, not highly oxygenated"
    
    # Look for methyl groups at specific positions
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 3:
        return False, "Less than 3 methyl groups found"
    
    # Count rotatable bonds as a rough measure of size
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a triterpenoid"
    
    # If all checks pass, it's likely a limonoid
    return True, "Contains a furan ring and a highly oxygenated triterpene skeleton with methyl groups at key positions"
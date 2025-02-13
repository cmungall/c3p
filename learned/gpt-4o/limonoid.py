"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids potentially containing a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Look for furan ring; limonoids often include a furan or furan-like structure
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if furan_pattern and not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan-like ring found"
    
    # Count oxygens; limonoids are typically highly oxygenated
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Insufficient oxygenation, found {o_count} oxygens"
    
    # Count carbons; triterpenoid-based count range typically
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 40:
        return False, f"Unusual triterpenoid carbon count, found {c_count} carbons"
    
    # Check for polycyclic framework
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, f"Too few rings, found {n_rings}"
    
    # Check for 4,4,8-trimethyl-17-furanylsteroid skeleton-like feature
    skeleton_pattern = Chem.MolFromSmarts("C1(C)CC[C@@H]2C3C(=O)OC4(C)C[C@H](OC3=1)C3C(OC4=O)C23C[C@@H]1(furanyl)")
    if skeleton_pattern and not mol.HasSubstructMatch(skeleton_pattern):
        return False, "Missing backbone structure resembling 4,4,8-trimethyl-17-furanylsteroid"
    
    return True, "Matches limonoid structural features (furan ring, high oxygenation, polycyclic skeleton)"
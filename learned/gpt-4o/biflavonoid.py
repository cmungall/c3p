"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer consisting of at least two flavonoid
    units linked together usually via C-C or C-O-C bonds.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define generic aromatic flavonoid unit pattern
    # This is a rough pattern capturing 2 aromatic rings with a keto group
    flavonoid_pattern = Chem.MolFromSmarts("c1ccccc1-c2cc(=O)cc(O)c2")
    
    # Look for at least two flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    if len(flavonoid_matches) < 2:
        return False, "Too few flavonoid units found"
    
    # Check for linkage between flavonoid units typically present in biflavonoids
    # Check for existence of direct linkage: C-C or C-O-C between flavonoid units
    linkage_pattern = Chem.MolFromSmarts("c1ccccc1-[C,O]-a2aaaaa2")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    if len(linkage_matches) < 2:
        return False, f"No or insufficient linkages between flavonoid units"
    
    # Check aromatic carbon count and presence of oxo groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbon atoms; likely not a biflavonoid"
    if o_count < 4:
        return False, "Too few oxygen atoms; likely not a biflavonoid"
    
    return True, "Contains multiple flavonoid units linked in a manner consistent with biflavonoids"
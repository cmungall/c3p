"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:39572 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid with a prototypical 4,4,8-trimethyl-17-furanylsteroid skeleton.

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
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"
    
    # Look for triterpenoid backbone (4 fused rings)
    backbone_pattern = Chem.MolFromSmarts("C1CCC2CCCC3CCCC(C1)C23")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Triterpenoid backbone not found"
    
    # Check for high oxygenation (at least 5 oxygens)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 5:
        return False, "Not enough oxygens for limonoid (need at least 5)"
    
    # Check for methyl groups at specific positions (4,4,8-trimethyl)
    methyl_pattern = Chem.MolFromSmarts("[C]1(CCC2CCCC3CCCC(C2)C3C)C")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing methyl groups at specific positions"
    
    # Check molecular weight - limonoids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for limonoid"
    
    # Count rings - limonoids typically have >=5 rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 5:
        return False, "Too few rings for limonoid (need at least 5)"
    
    return True, "Contains furan ring and triterpenoid backbone with high oxygenation and specific methyl groups"
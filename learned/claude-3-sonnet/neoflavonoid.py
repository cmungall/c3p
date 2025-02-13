"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: CHEBI:75768 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    Neoflavonoids have a 1-benzopyran core with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 1-benzopyran core
    # [#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1 represents benzene ring
    # [#8]1~[#6]~[#6] represents the pyran oxygen and adjacent carbons
    benzopyran_pattern = Chem.MolFromSmarts("[#8]1-[#6]~[#6]-[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2-1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran core found"

    # Look for aryl substituent at position 4
    # This pattern looks for an aromatic ring connected to the 4-position of the benzopyran
    aryl_4_pattern = Chem.MolFromSmarts("[#8]1-[#6]~[#6](-[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2)-[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2-1")
    if not mol.HasSubstructMatch(aryl_4_pattern):
        return False, "No aryl substituent at position 4"

    # Count number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Check for common substituents often found in neoflavonoids
    # Look for hydroxy, methoxy, or other oxygen-containing groups
    substituents_pattern = Chem.MolFromSmarts("([OH,OCH3,O])")
    substituent_matches = mol.GetSubstructMatches(substituents_pattern)
    
    # Most neoflavonoids have at least some oxygenated substituents
    if len(substituent_matches) == 0:
        return False, "No typical neoflavonoid substituents found"

    # Additional check for aromatic character
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "No aromatic systems found"

    # Check molecular weight - neoflavonoids typically >200 Da
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for neoflavonoid"

    return True, "Contains 1-benzopyran core with aryl substituent at position 4 and appropriate substituents"
"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a (potential) aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Must have exactly 5 carbons
    if c_count != 5:
        return False, f"Must have exactly 5 carbons, found {c_count}"
    
    # Should have 5 oxygens (4 hydroxyls + 1 ring oxygen or aldehyde oxygen)
    if o_count != 5:
        return False, f"Must have exactly 5 oxygens, found {o_count}"

    # Check for cyclic sugar patterns
    furanose_pattern = Chem.MolFromSmarts("[C]1[O][C]([C][C]1)")  # 5-membered ring
    pyranose_pattern = Chem.MolFromSmarts("[C]1[O][C]([C][C][C]1)")  # 6-membered ring
    
    # Check for aldehyde pattern (for open chain forms)
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    
    is_furanose = mol.HasSubstructMatch(furanose_pattern)
    is_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    if not (is_furanose or is_pyranose or has_aldehyde):
        return False, "Must be either a furanose, pyranose, or contain an aldehyde group"

    # Check for appropriate hydroxyl pattern
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Should have 4 or 5 hydroxyls (4 in cyclic form, 5 in open form)
    if hydroxyl_count not in [4, 5]:
        return False, f"Must have 4-5 hydroxyl groups, found {hydroxyl_count}"

    # Additional check for open chain form
    if has_aldehyde:
        # Verify it's a straight chain with hydroxyls
        chain_pattern = Chem.MolFromSmarts("[CH1](=O)[CH1][CH1][CH1][CH2]")
        if not mol.HasSubstructMatch(chain_pattern):
            return False, "Open chain form must have correct carbon backbone"
    
    # If we've made it here, it's an aldopentose
    form = "furanose" if is_furanose else "pyranose" if is_pyranose else "open chain"
    return True, f"Aldopentose in {form} form with correct number of carbons, oxygens, and hydroxyls"
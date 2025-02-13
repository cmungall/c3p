"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase attached to ribose/deoxyribose via N-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for furanose (5-membered sugar ring)
    furanose_pattern = Chem.MolFromSmarts("[CH2]-[CH1]-[CH1]-[CH1]-O1")
    if not mol.HasSubstructMatch(furanose_pattern):
        return False, "No ribose/deoxyribose sugar found"

    # Look for common nucleobase patterns
    base_patterns = [
        # Purine core
        "[nR1r6]1c[nR1r6]c2c1[nR1r6]c[nR1r6]c2",
        # Pyrimidine core 
        "[nR1r6]1c[nR1r6]cc(=O)[nR1r6]c1",
    ]
    
    found_base = False
    for pattern in base_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_base = True
            break
            
    if not found_base:
        return False, "No nucleobase pattern found"

    # Check for N-glycosidic bond (N-C1' bond between base and sugar)
    n_glycosidic = Chem.MolFromSmarts("[nR1r6]C1OC(CO)CC1")
    if not mol.HasSubstructMatch(n_glycosidic):
        return False, "No N-glycosidic bond found"

    # Basic size checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for nucleosides"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 9:
        return False, "Too few carbons for nucleoside"
    if n_count < 2:
        return False, "Too few nitrogens for nucleoside"
    if o_count < 4:
        return False, "Too few oxygens for nucleoside"

    # Check for reasonable ring count (typically 2-3 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Too few rings for nucleoside structure"
    if ring_info.NumRings() > 4:
        return False, "Too many rings for typical nucleoside"

    return True, "Contains nucleobase attached to ribose/deoxyribose via N-glycosidic bond"
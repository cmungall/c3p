"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and vary 
    only in carbon side chains and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid core (4 fused rings)
    steroid_core = Chem.MolFromSmarts("[CH2,CH]1[CH2,CH]2[CH2,CH][CH2,CH][C,c]3[CH2,CH][CH2,CH][C,c]4[C,c][C,c][C,c][C,c]4[C,c]3[C,c][C,c]2[CH2,CH]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for hydroxyl group (typically at C3 position)
    hydroxyl = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatches(hydroxyl):
        return False, "No hydroxyl group found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Phytosterols typically have 27-30 carbons
    if c_count < 27 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical range for phytosterols (27-35)"

    # Check molecular weight - phytosterols typically 400-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 600:
        return False, "Molecular weight outside typical range for phytosterols"

    # Check for aliphatic side chain
    aliphatic_chain = Chem.MolFromSmarts("[CH2,CH3][CH2,CH][CH2,CH][CH,C]")
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No suitable aliphatic side chain found"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Check for characteristics that would exclude it from being a phytosterol
    # Like presence of unusual elements
    allowed_elements = {1, 6, 8} # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Contains elements other than C, H, and O"

    # If all checks pass, it's likely a phytosterol
    return True, "Contains steroid core, hydroxyl group, and appropriate side chain characteristic of phytosterols"
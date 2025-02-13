"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True/False for classification, reason for the classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core structure pattern for dihydroagarofuran skeleton
    # Represents the characteristic tricyclic system with specific connectivity
    core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C][C][C]([C])[C]3(O[C]1([C])[C])[C]2")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroagarofuran core structure"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring count for dihydroagarofuran skeleton"

    # Check for ester groups (typically present in these compounds)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches < 2:
        return False, "Insufficient ester groups"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Should have at least 15 carbons (sesquiterpenoid core)
    if c_count < 15:
        return False, "Insufficient carbon count for sesquiterpenoid"
    
    # Should have multiple oxygens due to ester groups
    if o_count < 4:
        return False, "Insufficient oxygen count"

    # Check for characteristic bridged ring system
    bridged_pattern = Chem.MolFromSmarts("[C]12[C][C][C]1[C][C]2")
    if not mol.HasSubstructMatch(bridged_pattern):
        return False, "Missing characteristic bridged ring system"

    # Check for typical molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for dihydroagarofuran sesquiterpenoids"

    # Count sp3 carbons (should have several)
    sp3_c = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_c < 8:
        return False, "Insufficient sp3 carbons for dihydroagarofuran skeleton"

    return True, "Contains dihydroagarofuran skeleton with characteristic features"
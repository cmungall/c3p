"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: flavin compounds
Based on dimethylisoalloxazine skeleton with N10 substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavins are derivatives of dimethylisoalloxazine with N10 substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core tricyclic pattern with N10
    # Simplified but specific pattern for the key structural features
    core_pattern = Chem.MolFromSmarts('[#7]1c2c(nc(=O)[nH]c2=O)n([#6;!H3])c2cc([#6])c([#6])cc12')
    if core_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    if not mol.HasSubstructMatch(core_pattern):
        # Try alternate oxidation state pattern
        alt_pattern = Chem.MolFromSmarts('[#7]1c2c([nH]c(=O)[nH]c2=O)n([#6;!H3])c2cc([#6])c([#6])cc12')
        if alt_pattern is None or not mol.HasSubstructMatch(alt_pattern):
            return False, "Missing dimethylisoalloxazine core structure"

    # Count key features
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient number of nitrogen atoms (needs at least 4)"

    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Missing carbonyl groups"

    # Check for proper ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring system (needs 3 fused rings)"

    # Verify N10 substitution
    n10_pattern = Chem.MolFromSmarts('[#6]1:c:c([#6]):c([#6]):c:c1[#7]([#6;!H3])')
    if n10_pattern is None or not mol.HasSubstructMatch(n10_pattern):
        return False, "Missing or incorrect N10 substitution"

    # Additional check for aromatic system
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 6:
        return False, "Insufficient aromatic system"

    return True, "Contains dimethylisoalloxazine core with N10 substitution"

def test_flavins():
    """Test function with known flavins"""
    test_cases = [
        # Riboflavin
        "CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1",
        # FMN
        "C12=NC(NC(C1=NC=3C(N2C[C@@H]([C@@H]([C@@H](COP(=O)(O)O)O)O)O)=CC(=C(C3)C)C)=O)=O",
        # Lumiflavin
        "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C",
        # Not a flavin (benzene)
        "c1ccccc1",
        # Not a flavin (purine)
        "n1c2[nH]cnc2nc1"
    ]
    
    for smiles in test_cases:
        result, reason = is_flavin(smiles)
        print(f"SMILES: {smiles}")
        print(f"Is flavin: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_flavins()
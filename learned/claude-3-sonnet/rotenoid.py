"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: rotenoid compounds
Based on tetrahydrochromenochromene skeleton with specific substitution patterns
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure patterns
    # Basic chromeno-chromene skeleton
    core_pattern = Chem.MolFromSmarts("O1c2ccccc2C(=O)C2=C1CCOc1ccccc12")
    
    # Alternative core pattern with different saturation
    core_pattern2 = Chem.MolFromSmarts("O1c2ccccc2C3OC4=C(C3=O)c3ccccc3O4")
    
    # Check for presence of core structure
    if not (mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(core_pattern2)):
        return False, "Missing core chromenochromene skeleton"

    # Count oxygen atoms - rotenoids typically have multiple oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Too few oxygen atoms for rotenoid structure"

    # Look for carbonyl group (typically present in rotenoids)
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing carbonyl group characteristic of rotenoids"

    # Check for common substituents (methoxy groups are common)
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))

    # Count rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, "Insufficient ring count for rotenoid structure"

    # Additional structural features that might be present
    features = []
    if methoxy_count > 0:
        features.append(f"contains {methoxy_count} methoxy groups")
    
    # Check for specific ring patterns often found in rotenoids
    dioxole_pattern = Chem.MolFromSmarts("OCOCO")
    if mol.HasSubstructMatch(dioxole_pattern):
        features.append("contains dioxole ring")

    feature_str = " and ".join(features)
    if feature_str:
        feature_str = f", {feature_str}"

    return True, f"Contains chromenochromene core skeleton with required oxygenation pattern{feature_str}"
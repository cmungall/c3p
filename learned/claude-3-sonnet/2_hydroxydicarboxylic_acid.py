"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:35681 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid has two carboxylic acid groups and a hydroxy group
    on the alpha carbon of one of the carboxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxylic acid groups, need exactly 2"

    # Find alpha-hydroxy carboxylic acid pattern
    # Carbon with OH group that's also connected to a carboxyl group
    alpha_hydroxy_pattern = Chem.MolFromSmarts('[OH1]-[CH1,CH0]-C(=O)[OH1]')
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    
    if not alpha_hydroxy_matches:
        return False, "No alpha-hydroxy carboxylic acid group found"

    # Get total atom count (excluding H)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    
    # Check if molecule is too complex (likely not a simple dicarboxylic acid)
    if heavy_atom_count > 15:  # Most 2-hydroxydicarboxylic acids are relatively small
        return False, "Molecule too complex for a simple 2-hydroxydicarboxylic acid"

    # Count number of carbons and oxygens to ensure reasonable composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count != 5:  # 2 from carboxylic acids (2*2=4) + 1 from hydroxy = 5
        return False, "Must have exactly 5 oxygen atoms"

    # Check for presence of unexpected elements
    allowed_elements = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Contains elements other than C, H, and O"

    # Check ring count - most simple 2-hydroxydicarboxylic acids are not cyclic
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings - most 2-hydroxydicarboxylic acids are acyclic"

    return True, "Contains two carboxylic acid groups with a hydroxy group on the alpha carbon"
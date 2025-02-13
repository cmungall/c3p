"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:27102 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzoquinone core (cyclohexa-2,5-diene-1,4-dione)
    benzoquinone_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No benzoquinone core found"

    # Look for two methoxy groups (-OC) adjacent to quinone
    dimethoxy_pattern = Chem.MolFromSmarts("O=C1C(OC)=C(OC)C(=O)C=C1")
    if not mol.HasSubstructMatch(dimethoxy_pattern):
        return False, "Missing required 2,3-dimethoxy pattern"

    # Count methoxy groups (-OC) to ensure exactly 2
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_matches = len(mol.GetSubstructMatches(methoxy_pattern))
    if methoxy_matches < 2:
        return False, "Insufficient methoxy groups"

    # Look for at least one methyl group on the quinone ring
    methyl_quinone_pattern = Chem.MolFromSmarts("O=C1C(OC)=C(OC)C(=O)C(C)=C1")
    if not mol.HasSubstructMatch(methyl_quinone_pattern):
        return False, "Missing methyl group on quinone ring"

    # Optional: Check for isoprenoid/polyprenoid side chain
    # Pattern for isoprenoid unit CC=C(C)C or CC=C(C)CC
    isoprenoid_pattern = Chem.MolFromSmarts("CC=C(C)[C,CC]")
    has_isoprenoid = mol.HasSubstructMatch(isoprenoid_pattern)

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count < 4:
        return False, "Insufficient oxygen atoms for ubiquinone"
    
    # Basic structure confirmed
    if has_isoprenoid:
        return True, "Contains benzoquinone core with 2,3-dimethoxy groups, methyl group, and isoprenoid side chain"
    else:
        return True, "Contains benzoquinone core with 2,3-dimethoxy groups and methyl group"
"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:26611 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for uronic acid units (e.g., glucuronic acid, iduronic acid)
    uronic_acid_pattern = Chem.MolFromSmarts("[O-]C=C[C@H](O)[C@H](O)[C@H](O)CO")
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    
    # Look for glycosamine units (e.g., N-acetylglucosamine, N-acetylgalactosamine)
    glycosamine_pattern = Chem.MolFromSmarts("[NH2][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    
    # Check if uronic acids and glycosamines alternate
    if not uronic_acid_matches or not glycosamine_matches:
        return False, "Missing uronic acid or glycosamine units"
    
    # Look for sulfate groups
    sulfate_pattern = Chem.MolFromSmarts("[O-,OS(=O)(=O)O]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    
    # Count sulfur atoms to estimate degree of sulfation
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # Classify as mucopolysaccharide if uronic acids, glycosamines, and some sulfation
    if uronic_acid_matches and glycosamine_matches and (sulfate_matches or s_count > 0):
        return True, "Contains alternating uronic acid and glycosamine units, partially esterified with sulfate groups"
    else:
        return False, "Does not match mucopolysaccharide structural features"
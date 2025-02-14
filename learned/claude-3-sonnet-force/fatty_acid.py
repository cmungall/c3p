"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:36976 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid derived from or contained in esterified form
    in an animal or vegetable fat, oil or wax. Natural fatty acids commonly have a chain of 4 to 28
    carbons (usually unbranched and even-numbered), which may be saturated or unsaturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[C;$(C(=O)(O))]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for aliphatic chain (linear, branched, or cyclic)
    aliphatic_pattern = Chem.MolFromSmarts("[C;H3,H2,H1]~[C;H3,H2,H1]~[C;H3,H2,H1]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if not aliphatic_matches:
        return False, "No aliphatic chain found"
    
    # Check chain length (4 to 28 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4 or carbon_count > 28:
        return False, "Chain length outside the typical range for fatty acids"
    
    # Check for common functional groups
    functional_groups = ['[OX1]', '[OX2H]', '[OX2]', '[OX3H]', '[OX3]', '[Br]', '[Cl]', '[F]', '[I]', '[N]', '[S]', '[P]']
    functional_group_patterns = [Chem.MolFromSmarts(pattern) for pattern in functional_groups]
    has_functional_group = any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)
    
    # Check molecular weight (typically between 100 and 400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 400:
        return False, "Molecular weight outside the typical range for fatty acids"
    
    # Check for stereochemistry
    stereochem_pattern = Chem.MolFromSmarts("[C@@]")
    has_stereochem = mol.HasSubstructMatch(stereochem_pattern)
    
    # Classify as fatty acid
    reason = "Contains a carboxylic acid group and an aliphatic chain"
    if has_functional_group:
        reason += " with functional groups"
    if has_stereochem:
        reason += " and stereochemistry"
    
    return True, reason
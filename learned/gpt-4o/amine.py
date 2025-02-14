"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is derived from ammonia by replacing one, two, or three hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse Mol from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More refined patterns to detect amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary amine: NH2 bound to non-carbonyl
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1;!$(NC=O)]")   # Secondary amine: NH bound to non-carbonyl
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]")       # Tertiary amine: N bound to non-carbonyl
    cyclic_amine_pattern = Chem.MolFromSmarts("[NX3;R]")               # Cyclic amine: N in a ring

    # Check the molecule for primary, secondary, tertiary, or cyclic amine patterns
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Molecule contains a primary amine group"
    elif mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Molecule contains a secondary amine group"
    elif mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Molecule contains a tertiary amine group"
    elif mol.HasSubstructMatch(cyclic_amine_pattern):
        return True, "Molecule contains a cyclic amine group"

    return False, "No amine group detected"
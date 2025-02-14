"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:36285 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy (-OH) functional group
    in the alpha- or 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Look for hydroxy group (-OH) attached to C2
    hydroxy_pattern = Chem.MolFromSmarts("[CH2][CH]([OH])C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group at C2 position"

    # Check for long carbon chain (fatty acid)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long carbon chain found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 6:
        return False, "Too few carbons for a fatty acid"
    if o_count != 3:
        return False, "Must have exactly 3 oxygens (carboxyl and hydroxy groups)"

    return True, "Contains a long carbon chain with a carboxylic acid and hydroxy group at C2"
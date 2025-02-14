"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:18198 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclitol core pattern
    cyclitol_pattern = Chem.MolFromSmarts("[OX3r6,OX2r6][CX4r6][CX4r6][CX4r6][CX4r6][CX4r6]")
    if not mol.HasSubstructMatch(cyclitol_pattern):
        return False, "No cyclitol core found"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for correct atom counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if c_count != 7 or o_count != 7 or h_count != 14:
        return False, "Incorrect atom counts for quinic acid"

    # Look for common substituents (caffeoyl, feruloyl, coumaroyl)
    sub_patterns = [Chem.MolFromSmarts("O=C/C=C/c:1:c:c:c:c:c:1")]  # Caffeoyl
    sub_patterns.append(Chem.MolFromSmarts("COc:1:c:c:c:c:c:1/C=C/C=O"))  # Feruloyl
    sub_patterns.append(Chem.MolFromSmarts("O=C/C=C/c:1:c:c:c:c:c:1"))  # Coumaroyl

    sub_matches = [mol.GetSubstructMatches(p) for p in sub_patterns]
    num_subs = sum(len(matches) for matches in sub_matches)

    if num_subs > 3:
        return False, "Too many substituents for quinic acid"

    return True, "Contains cyclitol core with carboxylic acid group"
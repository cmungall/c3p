"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:36976 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) True if molecule is an amino sugar, False otherwise.
               Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sugar backbone
    sugar_pattern = Chem.MolFromSmarts("[OX2,OX1][CX4][CX4][CX4][CX4][OX2,OX1]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar backbone found"

    # Check for amino groups
    amino_pattern = Chem.MolFromSmarts("[NX3]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino groups found"

    # Check for hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No hydroxy groups found"

    # Additional checks
    n_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if n_oxygen < 3 or n_nitrogen < 1 or n_carbon < 3:
        return False, "Does not meet minimum atom count requirements"

    return True, "Contains a sugar backbone with at least one amino group"
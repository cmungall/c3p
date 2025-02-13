"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:33899 secondary amine
"""
from typing import Tuple
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def is_secondary_amine(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a compound formally derived from ammonia by replacing
    two hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atom found"

    # Check for secondary amine pattern
    sec_amine_pattern = Chem.MolFromSmarts("N($(*C),$(*C))")
    sec_amine_matches = mol.GetSubstructMatches(sec_amine_pattern)

    if len(sec_amine_matches) == 0:
        return False, "No secondary amine pattern found"

    # Check for other nitrogen patterns (primary, tertiary, quaternary)
    other_nitrogen_patterns = [
        Chem.MolFromSmarts("N($(*C),$(*N))"),  # Tertiary amine
        Chem.MolFromSmarts("N($(*C),$(*C),$(*C),$(*C))"),  # Quaternary amine
        Chem.MolFromSmarts("N($(*C))"),  # Primary amine
    ]

    for pattern in other_nitrogen_patterns:
        other_matches = mol.GetSubstructMatches(pattern)
        if len(other_matches) > 0:
            return False, "Molecule contains other nitrogen patterns"

    return True, "Contains secondary amine group(s)"
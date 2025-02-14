"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMART pattern for amines (includes primary, secondary, and tertiary setups)
    amine_pattern_primary_secondary = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')
    amine_pattern_tertiary = Chem.MolFromSmarts('[NX3;H0;!$(NC=O)]')

    # Count amino groups (primary and secondary amines, excluding amides)
    amine_matches_ps = mol.GetSubstructMatches(amine_pattern_primary_secondary)
    amine_matches_t = mol.GetSubstructMatches(amine_pattern_tertiary)

    amine_count_ps = len(amine_matches_ps)
    amine_count_t = len(amine_matches_t)

    # Total amine count considering primary, secondary, and tertiary nitrogens
    amine_count = amine_count_ps + amine_count_t

    # Determine classification
    if amine_count >= 2:
        return True, f"Molecule has {amine_count} amino groups, classifying as polyamine"
    else:
        return False, f"Molecule has {amine_count} amino group(s), not enough to classify as polyamine"
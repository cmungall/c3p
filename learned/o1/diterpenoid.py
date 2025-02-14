"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:35196 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is derived from a diterpene (composed of four isoprene units),
    but may have rearranged or modified skeletons, often missing some methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Extract counts of C, H, O, and other elements
    from collections import Counter
    import re

    # Use regular expressions to extract element counts from the formula
    elements = re.findall('([A-Z][a-z]*)(\d*)', formula)
    element_counts = Counter()
    for elem, count in elements:
        count = int(count) if count else 1
        element_counts[elem] += count

    c_count = element_counts.get('C', 0)
    h_count = element_counts.get('H', 0)
    o_count = element_counts.get('O', 0)
    n_count = element_counts.get('N', 0)
    s_count = element_counts.get('S', 0)
    other_elements = set(element_counts.keys()) - {'C', 'H', 'O', 'N', 'S'}

    # Diterpenoids generally have carbon counts around 20 (from four isoprene units)
    if c_count < 16 or c_count > 24:
        return False, f"Carbon count ({c_count}) not typical for diterpenoids (16-24 carbons)"

    # Calculate Double Bond Equivalents (DBE)
    dbe = (2 * c_count + 2 + n_count - h_count - halogen_count(mol)) / 2
    if dbe < 4:
        return False, f"DBE ({dbe}) too low for diterpenoids"

    # Check for presence of other elements
    if other_elements:
        return False, f"Contains elements not typical in diterpenoids: {', '.join(other_elements)}"

    # Check for terpenoid skeletal patterns (acyclic and cyclic diterpenes)
    # Common cyclic diterpene skeletons can be matched using SMARTS patterns
    # Here, we consider several core skeletons of diterpenes

    # List of common diterpene skeleton SMARTS patterns
    diterpene_skeletons = [
        # Labdane skeleton
        "[C&H1]1([C&H2])CC[C&H2]2[C&H2]1CC[C&H2]2",
        # Clerodane skeleton
        "[C&H1]1([C&H2])CC[C&H2]2[C&H2]1CC[C&H2]2C",
        # Abietane skeleton
        "C1=CC[C&H2]2[C&H2]1CCC3C2CCC=C3",
        # Gibberellane skeleton
        "[C&H2]1CC[C&H2]2[C&H2]1CC[C&H2]3[C&H2]2CC[C&H2]3",
        # Kaurane skeleton
        "[C&H2]1CCC2[C&H2]1CC[C&H2]3[C&H2]2CCC[C&H2]3",
        # Retinol skeleton (acyclic diterpene)
        "C=C(C)CCC=C(C)CCC=C(C)CC=C(C)C",
    ]

    matches_skeleton = False
    for pattern in diterpene_skeletons:
        skeleton = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(skeleton):
            matches_skeleton = True
            break

    if not matches_skeleton:
        return False, "Does not match common diterpene skeletons"

    # Check for rearranged or modified skeletons by looking for isoprene units
    # Count the number of C5 units
    c5_fragments = Chem.MolFromSmarts('[C&R]1~[C&R]~[C&R]~[C&R]~[C&R]1')
    c5_matches = mol.GetSubstructMatches(c5_fragments)
    if len(c5_matches) < 2:
        return False, f"Found {len(c5_matches)} C5 units, less than expected for diterpenoids"

    # Check if molecule is too oxidized or contains atypical functional groups
    # Diterpenoids may contain hydroxyl, carbonyl, carboxyl groups but rarely nitro groups, sulfates, etc.
    # List of atypical functional groups
    atypical_groups = [
        Chem.MolFromSmarts('[N+](=O)[O-]'),  # Nitro group
        Chem.MolFromSmarts('S(=O)(=O)[O-]'),  # Sulfate group
        Chem.MolFromSmarts('P(=O)(O)(O)O'),  # Phosphate group
    ]
    for group in atypical_groups:
        if mol.HasSubstructMatch(group):
            return False, "Contains atypical functional groups for diterpenoids"

    return True, "Molecule has features consistent with diterpenoids (carbon count, DBE, terpenoid skeleton)"

def halogen_count(mol):
    """
    Counts the number of halogen atoms in the molecule.

    Args:
        mol (rdkit.Chem.Mol): Molecule object

    Returns:
        int: Number of halogen atoms
    """
    halogens = [9, 17, 35, 53, 85]  # Atomic numbers for F, Cl, Br, I, At
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in halogens)
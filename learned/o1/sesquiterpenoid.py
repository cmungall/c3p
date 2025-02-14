"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is a terpenoid derived from a sesquiterpene, built from three isoprene units (C5 units),
    typically containing around 15 carbons, but may have rearrangements or modifications by the removal or addition of skeletal atoms.

    Due to the diversity and complexity of sesquiterpenoid structures, this function attempts to identify common sesquiterpene cores
    and patterns. It is not guaranteed to be accurate for all sesquiterpenoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Sesquiterpenoids typically have around 15 carbons, but can vary due to modifications
    if c_count < 12 or c_count > 25:
        return False, f"Carbon count is {c_count}, which is atypical but possible for a sesquiterpenoid"

    # Define common sesquiterpene skeletons (simplified)
    sesquiterpene_smarts = [
        # Farnesane skeleton (linear sesquiterpene)
        "C(C)(C)CC(C)(C)CC(C)(C)C",
        # Germacrane skeleton (cyclic sesquiterpene)
        "C1CCC(C=C1)CC=C(C)C",
        # Eudesmane skeleton
        "C1CC2CCCCC2(C)C1C",
        # Guaiane skeleton
        "C1CC2CC1CCC2C(C)C",
        # Cadinane skeleton
        "C1CC2C(C)CCC1CCC2",
        # Bisabolane skeleton
        "CC(C)CCCC(C)CC=C",
        # Other general sesquiterpene patterns
    ]

    # Check for matches to common sesquiterpene skeletons
    skeleton_found = False
    for smarts in sesquiterpene_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            skeleton_found = True
            break

    if not skeleton_found:
        # Sesquiterpenoids may have rearranged skeletons, so proceed to other checks
        pass

    # Check for terpenoid-like features: presence of isoprene units
    # An isoprene unit is C5H8, which can be challenging to identify due to rearrangements
    # Attempt to find three isoprene-like units by counting methyl (CH3) groups attached to a carbon chain
    methyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1)
    
    if methyl_count < 3:
        return False, f"Found {methyl_count} methyl groups, less than expected for a sesquiterpenoid"

    # Check for terpenoid functional groups (more specific)
    functional_groups = {
        'alpha-methylene-gamma-lactone': 'O=C1OC(CC=C1)=C',
        'alpha,beta-unsaturated ketone': '[C;!R]=[C;!R]-[C;!R]=O',
        'epoxide': 'C1OC1',  # Epoxide ring
        'furan': 'c1ccoc1',  # Furan ring
        'alcohol': '[CX4][OX2H]',  # Aliphatic alcohol
        'phenol': 'c[OH]',  # Phenolic hydroxyl
    }

    fg_matches = 0
    for fg_name, fg_smarts in functional_groups.items():
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if fg_pattern and mol.HasSubstructMatch(fg_pattern):
            fg_matches += 1

    if fg_matches == 0:
        return False, "No specific terpenoid functional groups found"

    # Allow heteroatoms since sesquiterpenoids may contain O, N, S, Cl, etc.
    # No need to exclude molecules with heteroatoms

    # If some criteria are met, classify as sesquiterpenoid
    if skeleton_found or fg_matches > 0:
        return True, "Molecule meets criteria for a sesquiterpenoid (skeleton match or functional groups)"

    return False, "Molecule does not meet criteria for a sesquiterpenoid"
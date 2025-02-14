"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""

from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove hydrogens for accurate atom counting
    mol = Chem.RemoveHs(mol)

    # Check for invalid atoms (exclude metals and halogens)
    allowed_atomic_nums = {1, 6, 7, 8}  # Allow C, H, N, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains invalid atom type {atom.GetSymbol()}"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Molecule has {c_count} carbon atoms, expected 6"

    # Exclude molecules with carboxylic acid groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Molecule contains carboxylic acid group, not a hexose"

    # Exclude molecules with ester groups (including lactones)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Molecule contains ester group (possible lactone), not a hexose"

    # Exclude molecules with phosphate or sulfate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    if mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(sulfate_pattern):
        return False, "Molecule contains phosphate or sulfate group, not a hexose"

    # Exclude molecules with large substituents (e.g., aromatic rings)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Molecule contains aromatic rings, not characteristic of hexoses"

    # Check for aldehyde group at position 1 (aldohexose)
    aldehyde_pattern = Chem.MolFromSmarts("[#6H1](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Check for ketone group at position 2 (ketohexose)
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # Check for cyclic hemiacetal (pyranose form) or hemiketal (furanose form)
    hemiacetal_pattern = Chem.MolFromSmarts("[C;R][O;R][C;R]")
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)

    # Ensure molecule is a monosaccharide (no glycosidic bonds)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;R1]-[O]-[C;!R]")
    if mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "Molecule appears to have a glycosidic bond, may not be a monosaccharide"

    # Evaluate if molecule has aldehyde, ketone, or cyclic form
    if has_aldehyde or has_ketone or has_hemiacetal:
        return True, "Molecule is a hexose: 6 carbons and contains aldehyde/ketone or cyclic hemiacetal/hemiketal"
    else:
        return False, "Molecule does not contain aldehyde, ketone, or cyclic hemiacetal/hemiketal in expected positions"
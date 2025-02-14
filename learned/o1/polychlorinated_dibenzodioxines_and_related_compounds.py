"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds

This program checks if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl, or related compound
by identifying core structures and counting chlorine or bromine substituents.

Improvements made:
- Updated SMARTS patterns for core structures using canonical SMILES to improve matching accuracy.
- Removed limits on the number of rings and heavy atoms to include larger structurally related molecules.
- Excluded molecules with hydroxyl groups by limiting allowed heteroatoms to oxygen in ring structures only.
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl,
    or related compound based on its SMILES string.

    These compounds are characterized by the presence of
    dibenzo-p-dioxin, dibenzofuran, or biphenyl core structures, and multiple chlorine or bromine atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is in the class, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure patterns using canonical SMILES
    # Dibenzo-p-dioxin core
    dibenzo_dioxin_smiles = 'O1c2ccccc2Oc3ccccc13'
    dibenzo_dioxin_pattern = Chem.MolFromSmiles(dibenzo_dioxin_smiles)

    # Dibenzofuran core
    dibenzofuran_smiles = 'O1c2ccccc2c3ccccc13'
    dibenzofuran_pattern = Chem.MolFromSmiles(dibenzofuran_smiles)

    # Biphenyl core
    biphenyl_smiles = 'c1ccccc1c2ccccc2'
    biphenyl_pattern = Chem.MolFromSmiles(biphenyl_smiles)

    # Check for core structures
    has_dibenzo_dioxin = mol.HasSubstructMatch(dibenzo_dioxin_pattern)
    has_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)

    if has_dibenzo_dioxin:
        core_structure = 'dibenzo-p-dioxin'
    elif has_dibenzofuran:
        core_structure = 'dibenzofuran'
    elif has_biphenyl:
        core_structure = 'biphenyl'
    else:
        return False, "Molecule does not contain dibenzo-p-dioxin, dibenzofuran, or biphenyl core structures"

    # Count number of chlorine and bromine atoms
    num_cl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    num_br = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35)
    total_halogens = num_cl + num_br

    if total_halogens < 2:
        return False, f"Molecule contains {total_halogens} halogen atom(s); requires at least 2 for classification"

    # Ensure that chlorine and bromine are the only halogens present
    num_f = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
    num_i = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53)
    if num_f > 0 or num_i > 0:
        return False, "Molecule contains fluorine or iodine atoms; only chlorine and bromine are allowed"

    # Check for prohibited functional groups (e.g., hydroxyl groups)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')  # Hydroxyl group
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Molecule contains hydroxyl groups; these are not allowed"

    # Ensure that the molecule does not contain other heteroatoms besides oxygen
    allowed_atomic_numbers = {6, 1, 8, 17, 35}  # Carbon, hydrogen, oxygen, chlorine, bromine
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_numbers:
            return False, f"Molecule contains atom {atom.GetSymbol()}; only C, H, O, Cl, and Br are allowed"

    return True, f"Molecule contains {core_structure} core and is polychlorinated/polybrominated with {total_halogens} halogen atoms"
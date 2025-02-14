"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: CHEBI:33441 metal atom
"""
from rdkit import Chem

# List of metal atomic numbers (based on the periodic table)
METAL_ATOMIC_NUMBERS = set([3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118])

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is an atom of an element that exhibits typical metallic properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule has only one atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"

    # Get atomic number of the atom
    atom = mol.GetAtoms()[0]
    atomic_num = atom.GetAtomicNum()

    # Check if atomic number corresponds to a metal
    if atomic_num in METAL_ATOMIC_NUMBERS:
        # Additional checks for false positives
        if atom.GetTotalNumHs() > 0 or atom.GetFormalCharge() != 0:
            return False, f"Atomic number {atomic_num} corresponds to a metal atom, but bonding pattern or charge suggests a non-metal"
        else:
            return True, f"Atomic number {atomic_num} corresponds to a metal atom"
    else:
        return False, f"Atomic number {atomic_num} does not correspond to a metal atom"
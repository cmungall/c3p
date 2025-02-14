"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHEBI:33340 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen atom based on its SMILES string.
    A chalcogen is any neutral p-block element belonging to group 16 of the periodic table:
    Oxygen (O), Sulfur (S), Selenium (Se), Tellurium (Te), Polonium (Po), Livermorium (Lv).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neutral chalcogen atom, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the list of atoms in the molecule
    atoms = mol.GetAtoms()

    # Check if the molecule consists of a single atom
    if len(atoms) != 1:
        return False, "Molecule does not consist of a single atom"

    # Get the atom
    atom = atoms[0]

    # Get the atomic number
    atomic_num = atom.GetAtomicNum()
    # Get the formal charge
    formal_charge = atom.GetFormalCharge()

    # List of atomic numbers of chalcogens
    chalcogen_atomic_numbers = [8, 16, 34, 52, 84, 116]

    # Check if the atomic number is in the list of chalcogens
    if atomic_num not in chalcogen_atomic_numbers:
        return False, f"Atom is not a chalcogen element (atomic number {atomic_num})"

    # Check if the atom is neutral
    if formal_charge != 0:
        return False, f"Atom has a formal charge of {formal_charge}, not neutral"

    return True, f"Atom is a neutral chalcogen element with atomic number {atomic_num}"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33340',
        'name': 'chalcogen',
        'definition': 'Any p-block element belonging to the group 16 family of the periodic table.'
    }
}
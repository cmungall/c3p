"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    Chalcogens are elements of group 16: O, S, Se, Te, Po.
    Handles isotopes of these atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains a single atom
    if mol.GetNumAtoms() != 1:
         return False, "Molecule is not a single atom."

    # Get the atom from the molecule
    atom = mol.GetAtomWithIdx(0)

    # Check the atom has no formal charge
    if atom.GetFormalCharge() != 0:
        return False, "Chalcogen must be neutral."

    # Get atomic number
    atomic_num = atom.GetAtomicNum()
    isotope = atom.GetIsotope()

    # Check if it's a chalcogen
    chalcogen_atomic_numbers = [8, 16, 34, 52, 84]  # O, S, Se, Te, Po
    if atomic_num in chalcogen_atomic_numbers:
        if isotope == 0 or (isotope > 0): #any isotope is OK
            return True, "Is a chalcogen atom"
        else:
            return False, f"Not a chalcogen atom. Isotope {isotope} is invalid"
    else:
        return False, "Not a chalcogen atom"
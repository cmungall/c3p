"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHEBI:33373 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    A chalcogen is any p-block element belonging to the group 16 family of the periodic table.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if there is only one atom in the molecule
    if mol.GetNumAtoms() != 1:
        return False, "Chalcogen must be a single atom"
    
    # Get the atomic number of the single atom
    atom = mol.GetAtoms()[0]
    atomic_number = atom.GetAtomicNum()
    
    # Check if the atomic number corresponds to a chalcogen element
    chalcogens = [8, 16, 34, 52, 84, 118]  # Atomic numbers of chalcogen elements
    if atomic_number in chalcogens:
        element_symbol = atom.GetSymbol()
        return True, f"The single atom with symbol '{element_symbol}' is a chalcogen element"
    else:
        return False, "The single atom is not a chalcogen element"
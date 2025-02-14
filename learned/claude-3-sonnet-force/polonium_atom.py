"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: CHEBI:33649 polonium atom
A radioactive metallic element discovered in 1898 by Marie Sklodowska Curie and named after her home country, Poland (Latin Polonia).
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polonium atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule has only one atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    # Get the atomic number of the atom
    atom = mol.GetAtoms()[0]
    atomic_number = atom.GetAtomicNum()
    
    # Check if the atomic number corresponds to polonium (84)
    if atomic_number == 84:
        return True, "Polonium atom"
    else:
        return False, f"Not a polonium atom (atomic number {atomic_number})"
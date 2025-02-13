"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: CHEBI:33369 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a SMILES string represents a polonium atom.
    
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
    
    # Check number of atoms
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    
    # Check if it's polonium
    if atom.GetSymbol() != "Po":
        return False, "Not a polonium atom"
    
    # Check if it has a mass number
    mass_num = atom.GetIsotope()
    if mass_num == 0:
        return False, "No isotope mass specified"
        
    # Check if mass number is in valid range for polonium
    # Based on examples, valid range appears to be 190-218
    if mass_num < 190 or mass_num > 218:
        return False, f"Invalid mass number {mass_num} for polonium (valid range: 190-218)"
    
    return True, f"Valid polonium isotope with mass number {mass_num}"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33369',
        'name': 'polonium atom',
        'definition': 'A radioactive metallic element discovered in 1898 by Marie Sklodowska Curie and named after her home country, Poland (Latin Polonia).',
        'parents': []
    }
}
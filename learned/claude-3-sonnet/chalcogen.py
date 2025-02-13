"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: chalcogen
A chalcogen is any p-block element belonging to the group 16 family of the periodic table.
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity is a chalcogen based on its SMILES string.
    Only identifies pure chalcogen atoms, not compounds or ions containing chalcogens.
    
    Args:
        smiles (str): SMILES string of the chemical entity
        
    Returns:
        tuple: (bool, str) indicating if it's a chalcogen and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if it's a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    
    # List of chalcogen atomic numbers
    chalcogens = {
        8: "oxygen",
        16: "sulfur",
        34: "selenium",
        52: "tellurium",
        84: "polonium"
    }
    
    atomic_num = atom.GetAtomicNum()
    
    # Check if it's a chalcogen
    if atomic_num not in chalcogens:
        return False, f"Not a chalcogen (atomic number {atomic_num})"
    
    # Check for formal charge - pure atoms should have no charge
    if atom.GetFormalCharge() != 0:
        return False, "Not a pure atom (has formal charge)"
        
    # Check for explicit hydrogens
    if atom.GetNumExplicitHs() > 0:
        return False, "Not a pure atom (has explicit hydrogens)"
        
    # Check for implicit hydrogens
    if atom.GetNumImplicitHs() > 0:
        return False, "Not a pure atom (has implicit hydrogens)"
        
    # Check for bonds
    if len(atom.GetBonds()) > 0:
        return False, "Not a pure atom (has bonds)"
    
    # Get isotope information if present
    isotope = atom.GetIsotope()
    element_name = chalcogens[atomic_num]
    
    if isotope:
        return True, f"Chalcogen: {element_name}-{isotope}"
    else:
        return True, f"Chalcogen: {element_name}"


__metadata__ = {
    'chemical_class': {
        'name': 'chalcogen',
        'definition': 'Any p-block element belonging to the group 16 family of the periodic table.',
    }
}
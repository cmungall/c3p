"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: metal atom
"""
from rdkit import Chem
import re

def is_metal_atom(smiles: str):
    """
    Determines if a SMILES string represents a metal atom.
    Must be a neutral, unbonded metal atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    
    # List of metal elements
    metal_elements = {
        'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',  # Alkali metals
        'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',  # Alkaline earth metals
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  # First-row transition metals
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',  # Second-row transition metals
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  # Third-row transition metals
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',  # Fourth-row transition metals
        'Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po',  # Post-transition metals
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',  # Lanthanides
        'Ho', 'Er', 'Tm', 'Yb', 'Lu',  # More lanthanides
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',  # Actinides
        'Es', 'Fm', 'Md', 'No', 'Lr'  # More actinides
    }
    
    # Check if SMILES is valid
    if not isinstance(smiles, str):
        return False, "Input must be a string"
    
    # Strict SMILES format for single atom: [<optional_mass>Element]
    # Only allows numbers (isotope) and element symbol inside brackets
    pattern = r'^\[(\d+)?([A-Z][a-z]?)\]$'
    match = re.match(pattern, smiles)
    if not match:
        return False, "Not a single atom SMILES (should be in format [Element] or [massElement])"
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check number of atoms
    if mol.GetNumAtoms() != 1:
        return False, "Must contain exactly one atom"
    
    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()
    
    # Check if it's a metal
    if symbol not in metal_elements:
        return False, f"{symbol} is not a metal element"
    
    # Check for formal charge - must be neutral
    if atom.GetFormalCharge() != 0:
        return False, f"{symbol} has a formal charge, not a neutral atom"
    
    # Check for implicit or explicit hydrogens
    if atom.GetNumImplicitHs() > 0 or atom.GetNumExplicitHs() > 0:
        return False, f"{symbol} has bonded hydrogens, not a free atom"
    
    # Check for any bonds
    if atom.GetDegree() > 0:
        return False, f"{symbol} has bonds to other atoms"
        
    return True, f"Single atom of metal element {symbol}"
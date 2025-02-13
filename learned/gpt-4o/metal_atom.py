"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a standalone metal atom based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a standalone metal atom, False otherwise.
        str: Reason for classification.
    """
    # Verified list of known metal symbols
    metal_symbols = [
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti",
        "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb",
        "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
        "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
        "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np",
        "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Fl",
        "Lv"
    ]

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that there's only a single atom in the molecule
    if mol.GetNumAtoms() != 1:
        return False, "Not a standalone atom"

    # Get the atom information
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()
    atomic_num = atom.GetAtomicNum()

    # Check for neutrality (no charge)
    if atom.GetFormalCharge() != 0:
        return False, f"{symbol} is a charged atom"

    # Ensure the standalone atom is an uncharged metal
    if symbol in metal_symbols and atomic_num > 0:
        return True, f"{symbol} is a standalone neutral metal atom"

    return False, f"{symbol} is not a standalone neutral metal atom"
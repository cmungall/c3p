"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:51868 gas molecular entity
Definition: Any main group molecular entity that is gaseous at standard temperature and
pressure (STP; 0degreeC and 100 kPa).
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get molecular weight
    mol_wt = Descriptors.MolWt(mol)
    
    # Check for elements known to be gases at STP
    gas_elements = [1, 2, 5, 6, 7, 8, 9, 10, 17, 18, 36, 54, 86]  # H, He, B, C, N, O, F, Ne, Cl, Ar, Kr, Xe, Rn
    element_nums = set([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    
    # Check for small gaseous organic compounds
    if all(elem_num in [1, 5, 6, 7, 8, 9, 17] for elem_num in element_nums):
        if mol_wt < 80 and len(mol.GetAtoms()) <= 8:
            return True, "Small organic compound, likely gaseous at STP"
    
    # Check for inorganic compounds known to be gases
    inorganic_gases = ["[O-][O+]=O", "O=C=O", "[C-]#[O+]", "ClCl", "FF", "Cl[H]", "I[H]", "[N-][N+]#[N]",
                       "BrBr", "BrF", "BrF3", "BrF5", "ClF", "ClF3", "ClF5", "ClO2", "F2O", "FNO2", "N2O4", "NOF",
                       "O3", "OF2", "SF6", "SO2F2", "SO3"]
    if Chem.MolToSmiles(mol) in inorganic_gases:
        return True, "Inorganic compound known to be gaseous at STP"
    
    # Check for noble gases and their isotopes
    noble_gases = [2, 10, 18, 36, 54, 86]  # He, Ne, Ar, Kr, Xe, Rn
    if len(element_nums) == 1 and list(element_nums)[0] in noble_gases:
        return True, "Noble gas or its isotope, known to be gaseous at STP"
    
    # Check for compounds with small molecular weight and containing only gas elements
    if element_nums.issubset(set(gas_elements)) and mol_wt < 100:
        # Check for specific functional groups or substructures
        if not Chem.MolFromSmarts("[OD1]"):  # Exclude compounds with hydroxyl groups
            return True, "Molecule contains only elements known to be gases at STP, and has low molecular weight"
    
    # If none of the above conditions are met, assume not a gas
    return False, "Not recognized as a gas molecular entity"
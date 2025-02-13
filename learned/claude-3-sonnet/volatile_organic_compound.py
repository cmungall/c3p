"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:51564 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import Tuple

def is_volatile_organic_compound(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is defined as an organic compound with an initial boiling point <= 250°C at 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule is organic (contains carbon)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
        return False, "Molecule does not contain carbon, not an organic compound"
    
    # Estimate boiling point
    boiling_point = Descriptors.BouchardeauBouchardat(mol)
    if boiling_point is None:
        return None, None  # Unable to estimate boiling point, cannot classify
    
    # Check if boiling point is <= 250°C
    if boiling_point <= 250:
        return True, f"Estimated boiling point of {boiling_point:.2f}°C is <= 250°C"
    else:
        return False, f"Estimated boiling point of {boiling_point:.2f}°C is > 250°C"
"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:35941 nitrohydrocarbon

A nitrohydrocarbon is a C-nitro compound that is a hydrocarbon in which one or more of the
hydrogens has been replaced by nitro groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule contains only C, H, N, O
    allowed_atoms = [6, 1, 7, 8]  # C, H, N, O
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if any(num not in allowed_atoms for num in atom_nums):
        return False, "Contains atoms other than C, H, N, O"
    
    # Check if molecule contains at least one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro groups present"
    
    # Check if molecule has a hydrocarbon skeleton
    hydrocarbon_pattern = Chem.MolFromSmarts("[#6;R]")  # Carbon in a ring or chain
    if not mol.HasSubstructMatch(hydrocarbon_pattern):
        return False, "No hydrocarbon skeleton present"
    
    return True, "Molecule is a nitrohydrocarbon"
"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:50729 organoiodine compound
An organoiodine compound is a compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for iodine atoms
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    if not iodine_atoms:
        return False, "No iodine atoms found"

    # Check for carbon-iodine bonds
    c_i_bonds = [bond for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in [atom.GetIdx() for atom in iodine_atoms] and bond.GetEndAtomIdx() in [atom.GetIdx() for atom in mol.GetAromaticAtoms()]]
    if not c_i_bonds:
        return False, "No carbon-iodine bonds found"

    return True, "Molecule contains at least one carbon-iodine bond"
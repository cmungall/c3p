"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound should contain at least one carbon-iodine bond.

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

    # Iterate over bonds and check for carbon-iodine bond
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check if there's a C-I bond
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 53) or \
           (atom1.GetAtomicNum() == 53 and atom2.GetAtomicNum() == 6):
            return True, "Contains carbon-iodine bond"

    return False, "No carbon-iodine bond found"
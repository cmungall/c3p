"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine bond.

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
    
    # Iterate through bonds to find a carbon-iodine bond
    for bond in mol.GetBonds():
        # Get the beginning and ending atoms of the bond
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()

        # Check for a carbon-iodine bond (atomic numbers: 6 for carbon, 53 for iodine)
        if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 53) or \
           (begin_atom.GetAtomicNum() == 53 and end_atom.GetAtomicNum() == 6):
            return True, "Contains at least one carbon-iodine bond"
    
    return False, "No carbon-iodine bond found"
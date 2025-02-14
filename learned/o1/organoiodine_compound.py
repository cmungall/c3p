"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: organoiodine compound
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine (C-I) bond where the carbon is part of an organic framework.

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

    # Exclude molecules with less than 3 carbon atoms to avoid small molecules like iodoform
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 3:
        return False, f"Only {num_carbons} carbon atom(s); insufficient for organoiodine compound"

    # Find all carbon-iodine bonds where carbon is connected to at least one other carbon
    is_organoiodine = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check for C-I bond
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 53) or \
           (atom1.GetAtomicNum() == 53 and atom2.GetAtomicNum() == 6):
            carbon_atom = atom1 if atom1.GetAtomicNum() == 6 else atom2
            # Check if carbon is connected to at least one other carbon (excluding the iodine)
            neighbor_carbons = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom2.GetIdx()]
            if neighbor_carbons:
                is_organoiodine = True
                break  # Found a valid C-I bond
    if is_organoiodine:
        return True, "Contains C-I bond with carbon in organic framework"
    else:
        return False, "No suitable carbon-iodine bonds found"
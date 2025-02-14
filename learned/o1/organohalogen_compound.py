"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound
Definition: 'A compound containing at least one carbon-halogen bond (where X is a halogen atom).'
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond (C-X), where X is F, Cl, Br, or I.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define atomic numbers for carbon and halogens
    carbon_atomic_num = 6
    halogen_atomic_nums = [9, 17, 35, 53]  # F, Cl, Br, I

    # Check for carbon-halogen bonds
    has_carbon_halogen_bond = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atomic_nums = (atom1.GetAtomicNum(), atom2.GetAtomicNum())
        if (carbon_atomic_num in atomic_nums) and (any(h in atomic_nums for h in halogen_atomic_nums)):
            has_carbon_halogen_bond = True
            break  # No need to check further

    if has_carbon_halogen_bond:
        return True, "Contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bonds found"
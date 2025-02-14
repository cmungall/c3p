"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: CHEBI:35486 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound is a compound containing at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon-chlorine bonds
    has_c_cl_bond = False
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 17:
            has_c_cl_bond = True
            break
        elif begin_atom.GetAtomicNum() == 17 and end_atom.GetAtomicNum() == 6:
            has_c_cl_bond = True
            break

    if has_c_cl_bond:
        return True, "Contains at least one carbon-chlorine bond"
    else:
        return False, "Does not contain any carbon-chlorine bonds"
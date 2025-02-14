"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one carbon-chlorine bond.

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

    # Iterate over atoms to find chlorine atoms bonded to carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # 17 is the atomic number of chlorine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # 6 is the atomic number of carbon
                    return True, "Contains at least one carbon-chlorine bond"

    return False, "No carbon-chlorine bond found"
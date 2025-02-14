"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:33842 primary amine
A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get atoms with atomic number 7 (nitrogen)
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if there is at least one nitrogen atom
    if not n_atoms:
        return False, "No nitrogen atom found"

    # Check if any nitrogen atom has exactly 2 hydrogen atoms attached
    for n_atom in n_atoms:
        if sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1) == 2:
            # This is a primary amine nitrogen
            break
    else:
        return False, "No primary amine group found"

    # Check if there is at least one carbon atom bonded to the nitrogen
    for n_atom in n_atoms:
        if any(neighbor.GetAtomicNum() == 6 for neighbor in n_atom.GetNeighbors()):
            return True, "Contains a primary amine group (-NH2) bonded to a carbon"

    return False, "No carbon atom bonded to primary amine group"
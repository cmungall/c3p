"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate has the general formula R-N=C=S.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define isothiocyanate SMARTS pattern
    # Nitrogen double bonded to carbon, which is double bonded to sulfur
    isothiocyanate_pattern = Chem.MolFromSmarts("[NX2]=C=[SX1]")

    # Check for isothiocyanate group
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    if not matches:
        return False, "Does not contain isothiocyanate functional group (N=C=S)"

    # Ensure the nitrogen is attached to an R group (any atom except C or H in the N=C=S group)
    for match in matches:
        n_idx = match[0]  # Index of the nitrogen atom in the match
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Get neighbors of the nitrogen atom excluding the carbon in N=C=S
        neighbors = [atom for atom in n_atom.GetNeighbors() if atom.GetIdx() != match[1]]
        if neighbors:
            return True, "Contains isothiocyanate functional group (R-N=C=S)"
    return False, "Isothiocyanate group not attached to any R group"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:47905',
        'name': 'isothiocyanate',
        'definition': 'An organosulfur compound with the general formula R-N=C=S.'
    }
}
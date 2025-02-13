"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:33104 amine
A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if at least one nitrogen atom is present
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"

    # Check if any nitrogen atom is trivalent and not part of a ring or specific functional group
    for n_atom in nitrogen_atoms:
        # Skip nitrogen atoms that are part of a ring
        if n_atom.IsInRing():
            continue

        # Skip nitrogen atoms that are part of nitro groups
        if sum(1 for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() == 8) > 1:
            continue

        # Skip nitrogen atoms that are part of amides
        if sum(1 for nbr in n_atom.GetNeighbors() if (nbr.GetAtomicNum() == 8 and
                                                      sum(1 for nbr_nbr in nbr.GetNeighbors()
                                                          if nbr_nbr.GetAtomicNum() == 8) == 1)) > 0:
            continue

        # Check if trivalent and has at least one N-C bond
        if len(n_atom.GetNeighbors()) == 3 and any(nbr.GetAtomicNum() == 6 for nbr in n_atom.GetNeighbors()):
            return True, "Contains a trivalent nitrogen atom with at least one N-C bond, not part of a ring, nitro group, or amide"

    return False, "No nitrogen atom meets the criteria for an amine"
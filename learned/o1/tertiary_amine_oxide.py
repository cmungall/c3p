"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide has a nitrogen atom with a positive charge (+1),
    bonded to an oxygen atom with a negative charge (-1) and three organic groups
    (typically carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms to find nitrogen atoms with formal charge +1
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == +1:
            neighbors = atom.GetNeighbors()
            # Check if nitrogen is bonded to exactly four atoms
            if len(neighbors) != 4:
                continue
            # Find oxygen neighbors with formal charge -1
            oxygen_neighbors = [nbr for nbr in neighbors
                                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1]
            # There must be exactly one such oxygen neighbor
            if len(oxygen_neighbors) != 1:
                continue
            # The other three neighbors should be carbon atoms (organic groups)
            other_neighbors = [nbr for nbr in neighbors if nbr.GetIdx() != oxygen_neighbors[0].GetIdx()]
            if all(nbr.GetAtomicNum() == 6 for nbr in other_neighbors):
                return True, "Contains tertiary amine oxide functional group"
    return False, "Does not contain tertiary amine oxide functional group"

__metadata__ = {   'chemical_class': {   'name': 'tertiary amine oxide',
                                         'definition': 'An N-oxide where there are three organic groups bonded to the nitrogen atom.'}}
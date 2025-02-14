"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:48817 1,2,4-triazines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_4_triazine(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a compound with a 1,2,4-triazine skeleton, where nitrogen atoms
    replace carbon at positions 1, 2, and 4 of the core benzene ring structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 1,2,4-triazine core (6-membered ring with 3 nitrogen atoms)
    triazine_pattern = Chem.MolFromSmarts("c1ncncn1")
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "No 1,2,4-triazine core found"

    # Check the positions of nitrogen atoms
    atoms = mol.GetAtoms()
    n_positions = []
    for atom in atoms:
        if atom.GetAtomicNum() == 7:  # Nitrogen
            n_positions.append(atom.GetIdx())

    # Nitrogen atoms should be at positions 1, 2, and 4 of the core ring
    if sorted(n_positions) != [0, 1, 3]:
        return False, "Nitrogen atoms not at positions 1, 2, and 4 of the core ring"

    return True, "Contains a 1,2,4-triazine skeleton with nitrogen atoms replacing carbon at positions 1, 2, and 4 of the core ring"
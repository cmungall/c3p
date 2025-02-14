"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a six-membered aromatic ring with nitrogen atoms at positions 1, 2, and 4.

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

    # Define the 1,2,4-triazine substructure
    triazine_smiles = 'c1cnncn1'
    triazine_mol = Chem.MolFromSmiles(triazine_smiles)
    if triazine_mol is None:
        return False, "Error in defining 1,2,4-triazine substructure"

    # Check if the molecule contains the 1,2,4-triazine ring
    if mol.HasSubstructMatch(triazine_mol):
        return True, "Contains 1,2,4-triazine ring"
    else:
        return False, "Does not contain 1,2,4-triazine ring"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '1,2,4-triazines',
        'definition': 'Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon at positions 1, 2 and 4 of the core benzene ring structure.',
        'parents': []
    }
}
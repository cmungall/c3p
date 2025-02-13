"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:36357 polyamine
Any organic amino compound that contains two or more amino groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of amino groups
    amino_pattern = Chem.MolFromSmarts("[NH2]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    n_amino_groups = len(amino_matches)

    if n_amino_groups < 2:
        return False, "Must have at least 2 amino groups to be a polyamine"

    # Check if the molecule is organic (contains carbon)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count == 0:
        return False, "Not an organic compound (no carbon atoms)"

    # Check for common inorganic polyamines like hydrazine
    inorganic_pattern = Chem.MolFromSmarts("[!#6;!#7]")
    if mol.HasSubstructMatch(inorganic_pattern):
        return False, "Contains non-organic atoms, likely an inorganic polyamine"

    return True, f"Organic compound with {n_amino_groups} amino groups"
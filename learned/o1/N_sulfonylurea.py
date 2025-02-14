"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea is a urea in which one of the hydrogens attached to a nitrogen
    of the urea group is replaced by a sulfonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-sulfonylurea
    # N-sulfonylurea: O=C(N)N-S(=O)(=O)R
    n_sulfonylurea_pattern = Chem.MolFromSmarts('NC(=O)N[S](=O)(=O)[#6]')
    if n_sulfonylurea_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for N-sulfonylurea substructure
    if mol.HasSubstructMatch(n_sulfonylurea_pattern):
        return True, "Contains N-sulfonylurea moiety"
    else:
        return False, "Does not contain N-sulfonylurea moiety"

__metadata__ = {
    'chemical_class': {
        'name': 'N-sulfonylurea',
        'definition': 'A urea in which one of the hydrogens attached to a nitrogen of the urea group is replaced by a sulfonyl group. The N-sulfonylurea moiety is a key group in various herbicides, as well as in a number of antidiabetic drugs used in the management of type 2 diabetes mellitus.',
    }
}
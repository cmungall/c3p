"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: CHEBI:50123 N-sulfonylurea
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Look for urea group (-N(C=O)N-)
    urea_pattern = Chem.MolFromSmarts("[NH1]C(=O)N")
    if not mol.HasSubstructMatch(urea_pattern):
        return False, "No urea group found"

    # Look for sulfonyl group (-SO2-)
    sulfonylurea_pattern = Chem.MolFromSmarts("S(=O)(=O)")
    if not mol.HasSubstructMatch(sulfonylurea_pattern):
        return False, "No sulfonyl group found"

    # Check if the molecule contains an N-sulfonylurea group
    nsulfonylurea_pattern = Chem.MolFromSmarts("[NH1]C(=O)N~[S](=O)(=O)")
    if mol.HasSubstructMatch(nsulfonylurea_pattern):
        return True, "Contains an N-sulfonylurea group (urea with a sulfonyl group attached to a nitrogen)"

    return False, "No N-sulfonylurea group found"
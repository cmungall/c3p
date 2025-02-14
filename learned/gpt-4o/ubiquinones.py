"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone has a 2,3-dimethoxy-5-methylbenzoquinone core with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for benzoquinone core with 2,3-dimethoxy and 5-methyl groups
    benzoquinone_pattern = Chem.MolFromSmarts("COC1=C(C)C(=O)C=CC1=O")
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No 2,3-dimethoxy-5-methylbenzoquinone core found"
    
    # Check for long polyprenoid chain at position 6 (using multiple isoprene units)
    polyprenoid_pattern = Chem.MolFromSmarts("C=C(C)CC")
    polyprenoid_matches = mol.GetSubstructMatches(polyprenoid_pattern)
    if len(polyprenoid_matches) < 1:
        return False, "No polyprenoid side chain found"

    return True, "Contains ubiquinone core with polyprenoid side chain"
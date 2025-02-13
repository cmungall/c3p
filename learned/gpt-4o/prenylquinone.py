"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a basic quinone pattern (1,4-benzoquinone)
    quinone_pattern = Chem.MolFromSmarts("C1=CC(=O)C=CC1=O")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone moiety found"

    # Define a basic pattern for prenyl side-chains (isoprene units: C=C-C=C-C)
    prenyl_sidechain_pattern = Chem.MolFromSmarts("C=C-C-C=C")
    prenyl_matches = mol.GetSubstructMatches(prenyl_sidechain_pattern)
    if not prenyl_matches:
        return False, "No polyprenyl side-chain found"

    return True, "Contains a quinone moiety with a prenyl side-chain"
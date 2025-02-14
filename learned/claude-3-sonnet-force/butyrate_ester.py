"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: CHEBI:34835 butyrate ester
Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for butyrate ester substructure (butanoate group attached to oxygen)
    butanoate_pattern = Chem.MolFromSmarts("CCCC(=O)O*")
    is_butyrate_ester = mol.HasSubstructMatch(butanoate_pattern)

    if is_butyrate_ester:
        return True, "Contains a butanoate group attached to an oxygen atom"
    else:
        return False, "Does not contain a butanoate group attached to an oxygen atom"
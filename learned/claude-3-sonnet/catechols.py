"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:27027 catechol
Any compound containing an o-diphenol component
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechol(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component (a benzene ring
    with two hydroxyl groups attached to adjacent carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define o-diphenol substructure pattern
    o_diphenol_pattern = Chem.MolFromSmarts("c1c(O)c(O)ccc1")

    # Check if molecule contains the o-diphenol substructure
    if mol.HasSubstructMatch(o_diphenol_pattern):
        return True, "Contains an o-diphenol component"
    else:
        return False, "Does not contain an o-diphenol component"
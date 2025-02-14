"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:38181 catechol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechol(smiles: str):
    """
    Determines if a molecule contains an o-diphenol (catechol) component based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a catechol component, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol pattern
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cc1")
    
    # Check if the molecule contains the catechol pattern
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains an o-diphenol (catechol) component"
    else:
        return False, "Does not contain an o-diphenol (catechol) component"
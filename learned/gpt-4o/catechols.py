"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol contains an o-diphenol component (two hydroxyl groups on adjacent carbons on an aromatic ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broader SMARTS patterns for catechol moiety
    # Detecting o-dihydroxy groups (ortho-diphenol) in potentially various aromatic rings
    catechol_patterns = [
        Chem.MolFromSmarts("c1c(O)c(O)c[c,n]c1"),  # Basic ortho-diphenol pattern, allowing variable connectivity
        Chem.MolFromSmarts("c1(O)c(O)cc[n,c]c1"),  # Including nitrogen-containing heterocycles
        Chem.MolFromSmarts("c1(O)c(O)c[c,n]cc1"),  # Extended aromatic rings with ortho-dihydroxy
        Chem.MolFromSmarts("Oc1ccc(O)c[c,H]c1"),  # Generalized aromatic with ortho hydroxys
        Chem.MolFromSmarts("Oc1ccc2cc(O)c(O)cc2c1"),  # Biphenyl systems
    ]

    # Check each pattern
    for pattern in catechol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains o-diphenol component in an aromatic ring system"
    
    return False, "Does not contain o-diphenol component in an aromatic ring system"
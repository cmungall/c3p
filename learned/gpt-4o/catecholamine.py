"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine typically has a catechol structure with an ethylamine group (2-carbon chain leading to an amine)
    and may include various substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exact catechol pattern (benzene 1,2-diol)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol group (benzene-1,2-diol) found"
    
    # Look for ethylamine group with flexibility for substitution
    ethylamine_patterns = [
        Chem.MolFromSmarts("c1cc(O)c(O)cc1CCN"),
        Chem.MolFromSmarts("c1cc(O)c(O)cc1C(C)N"),
        Chem.MolFromSmarts("c1cc(O)c(O)cc1CC[NH2]"),
        Chem.MolFromSmarts("c1c(O)cc(O)cc1CCN([CX4])"),
        # More patterns can be added to match the variations in SMILES samples
    ]

    for pattern in ethylamine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains catechol group with an ethylamine linkage or variant"

    return False, "No ethylamine group connected to catechol found"
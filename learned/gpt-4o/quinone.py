"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    Quinones are characterized by a fully conjugated cyclic dione structure,
    derived from aromatic compounds by converting an even number of -CH= groups
    into -C(=O)- groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns for various fully conjugated cyclic dione structures
    quinone_patterns = [
        Chem.MolFromSmarts("c1c(=O)c([OH1])c([OH1])c(=O)c1"),  # Simple benzoquinone
        Chem.MolFromSmarts("c1cc(=O)c2c(c1)ccc(=O)c2"),         # Naphthoquinone-like
        Chem.MolFromSmarts("c1cc2c(c1)c(=O)c3ccccc3c2=O"),     # Anthraquinone-like
        Chem.MolFromSmarts("Oc1ccc2C(=O)C=CC(=O)c2c1"),        # Extended conjugation
        Chem.MolFromSmarts("[OH1]c1cc2C(=O)C=C(C1)C(=O)c3cc(C2)c4ccccc4c3=O")  # Complex cyclic with layers
    ]
    
    # Check for quinone patterns
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains fully conjugated cyclic dione structure typical of quinones"
    
    return False, "Does not contain the structural features typical of quinones"
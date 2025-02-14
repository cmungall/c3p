"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule contains an azole ring based on its SMILES string.
    An azole is a five-membered heteroarene containing nitrogen and possibly other
    heteroatoms like sulfur or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    azole_patterns = [
        Chem.MolFromSmarts("n1cccn1"),  # Pyrrole-like azole
        Chem.MolFromSmarts("n1c[nH]cc1"),  # Imidazole
        Chem.MolFromSmarts("n1ncc[nH]1"),  # 1,2,3-Triazole
        Chem.MolFromSmarts("n1ncnn1"),  # 1,2,4-Triazole
        Chem.MolFromSmarts("n1nncn1"),  # Tetrazole
        Chem.MolFromSmarts("c1ncn[nH]1"),  # Variant for 1,2,4-Triazole
        Chem.MolFromSmarts("c1nocn1"),  # Oxazole
        Chem.MolFromSmarts("c1nosc1"),  # Thiazole
        Chem.MolFromSmarts("c1nscn1")  # Another Thiazole variant
    ]

    # Check if any azole pattern is present in the molecule
    for pattern in azole_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains an azole ring based on defined patterns"

    return False, "No azole ring found"

# Example test (remove in production code)
# is_azole("CNC1=NC(=O)[C@@H](O1)[C@H](C)c1c[nH]c2ccccc12")  # Test with a known azole
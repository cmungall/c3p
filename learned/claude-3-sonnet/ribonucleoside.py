"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:26199 ribonucleoside
A ribonucleoside is any nucleoside where the sugar component is D-ribose.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ribose sugar
    ribose_pattern = Chem.MolFromSmarts("[CR1]([CR]([CR]([CR]([CR]([OR])[OR])[OR])[OR])[OR])")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar found"

    # Check for nucleobase
    nucleobase_patterns = [
        Chem.MolFromSmarts("[nr3r5]"),  # Pyrimidine
        Chem.MolFromSmarts("[nr5r5]"),  # Imidazole
        Chem.MolFromSmarts("[nr5r6]"),  # Purine
    ]
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Check for glycosidic bond between ribose and nucleobase
    ribonucleoside_pattern = Chem.MolFromSmarts("[CR1]([CR]([CR]([CR]([CR]([OR])[OR])[OR])[OR])[OR])[nr3r5,nr5r5,nr5r6]")
    if not mol.HasSubstructMatch(ribonucleoside_pattern):
        return False, "No glycosidic bond between ribose and nucleobase"

    return True, "Contains a ribose sugar and a nucleobase connected by a glycosidic bond"
"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone contains a naphthoquinone core with one or more hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for naphthoquinone core and hydroxyl group
    naphthoquinone_pattern = Chem.MolFromSmarts("C1=CC=CC2=C1C(=O)C=CC2=O")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    
    # Check for naphthoquinone core
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core found"
    
    # Check for at least one hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl group found"
    
    return True, "Contains naphthoquinone core with at least one hydroxyl group"

# Example usage:
smiles_example = "OC1=CC(=O)c2ccccc2C1=O"  # Lawsone
result, reason = is_hydroxynaphthoquinone(smiles_example)
print(f"Is hydroxynaphthoquinone: {result}, Reason: {reason}")
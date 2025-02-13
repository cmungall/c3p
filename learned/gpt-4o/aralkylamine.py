"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an aralkylamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for an aromatic group
    # Ensure the molecule contains an aromatic ring
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return False, "No aromatic group found"

    # Ensure the molecule contains an amine group (excluding quaternary amines)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;!$([N+])]")):
        return False, "No amine group found"

    # Check for aralkyl-amine linkages, using SMARTS to ensure linkage specifically involves
    # an aromatic group connected through an alkyl chain to the amine group
    aralkylamine_pattern = Chem.MolFromSmarts("[a][CX4][NX3;!$([N+])]")
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "Aromatic group not properly attached through an alkyl to amine"

    return True, "Contains an alkylamine with aromatic substitution"

# Testing the function with examples should follow to ensure it meets the criteria correctly.
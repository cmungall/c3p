"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is defined as any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS pattern
    # This pattern represents the flavanone skeleton without substitutions
    flavanone_smarts = 'O=C1CC(Oc2ccccc2)Cc3ccccc13'
    flavanone_pattern = Chem.MolFromSmarts(flavanone_smarts)
    if flavanone_pattern is None:
        return False, "Invalid flavanone SMARTS pattern"

    # Check if the molecule contains the flavanone core
    flavanone_matches = mol.GetSubstructMatches(flavanone_pattern)
    if not flavanone_matches:
        return False, "Molecule does not contain the flavanone skeleton"

    # Define the dihydroflavonol SMARTS pattern
    # This pattern includes a hydroxy group at position 3 of the heterocyclic ring
    dihydroflavonol_smarts = 'O[C@H]1C[C@@H](O)C(=O)c2ccccc12'
    dihydroflavonol_pattern = Chem.MolFromSmarts(dihydroflavonol_smarts)
    if dihydroflavonol_pattern is None:
        return False, "Invalid dihydroflavonol SMARTS pattern"

    # Check if the molecule contains the dihydroflavonol pattern
    dihydroflavonol_matches = mol.GetSubstructMatches(dihydroflavonol_pattern)
    if not dihydroflavonol_matches:
        return False, "Molecule does not have a hydroxy group at position 3 of the heterocyclic ring"

    return True, "Molecule is a dihydroflavonol"
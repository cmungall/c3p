"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic compound where a single carbon atom is attached to both
    an amino group and a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define the SMARTS pattern for a carbon attached to both hydroxyl and amino groups
    # [CX4;!$(C(=O))]: sp3-hybridized carbon not double-bonded to oxygen (exclude carbonyl carbons)
    # ([OH]): hydroxyl group attached to the carbon
    # ([NX3;H2,H1,H0]): trivalent nitrogen (amine group) with 0-2 hydrogens
    hemiaminal_pattern = Chem.MolFromSmarts('[CX4;!$(C(=O))]([OH])[NX3;H2,H1,H0]')

    # Search for substructure matches in the molecule
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    if matches:
        return True, "Contains a carbon atom attached to both hydroxyl and amino groups (hemiaminal)."
    else:
        return False, "Does not contain a hemiaminal functional group."
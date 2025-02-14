"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:29790 catechin

Catechins are members of the class of hydroxyflavan that have a flavan-3-ol skeleton
and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the catechin core scaffold
    catechin_core = Chem.MolFromSmarts("[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1ccc(O)c(O)c1")

    # Check if the molecule contains the catechin core
    if not mol.HasSubstructMatch(catechin_core):
        return False, "Molecule does not contain the catechin core scaffold"

    # Define other common catechin substructures
    galloyl_group = Chem.MolFromSmarts("C(=O)c1cc(O)c(O)c(O)c1")
    ester_group = Chem.MolFromSmarts("C(=O)O")
    methoxy_group = Chem.MolFromSmarts("OC")
    hydroxyl_group = Chem.MolFromSmarts("O")

    # Count the occurrences of these substructures
    n_galloyl = len(mol.GetSubstructMatches(galloyl_group))
    n_ester = len(mol.GetSubstructMatches(ester_group))
    n_methoxy = len(mol.GetSubstructMatches(methoxy_group))
    n_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_group))

    # Catechins typically have multiple hydroxyl groups and may have galloyl, ester, or methoxy substituents
    if n_hydroxyl < 3:
        return False, "Too few hydroxyl groups for a catechin"

    # Check for other common catechin features
    fmcs = rdFMCS.FindMCS([mol, catechin_core], atomsToUse=[0], completeRingsOnly=True)
    if fmcs.numBonds < 13:
        return False, "Insufficient ring system similarity to catechins"

    return True, "Molecule contains the catechin core scaffold and has typical catechin substituents"
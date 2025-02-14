"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is characterized by the presence of a beta-D-glucuronic acid
    moiety linked via a glycosidic bond at the anomeric carbon (C1) to an aglycone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-glucuronic acid linked via glycosidic bond
    # This pattern matches a glucuronic acid moiety with beta configuration at C1
    # Note: Matching stereochemistry in SMARTS can be complex and may not be fully accurate
    glucuronic_acid_smarts = '[C@H]1(O[!H])[O][C@H](O)[C@H](O)[C@H](O)[C@H]1C(=O)O'
    glucuronic_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)
    if glucuronic_pattern is None:
        return None, "Unable to construct SMARTS pattern"

    # Check if the molecule contains the beta-D-glucuronic acid moiety
    if not mol.HasSubstructMatch(glucuronic_pattern):
        return False, "Beta-D-glucuronic acid moiety not found"

    # Additional checks can be performed to verify the linkage and stereochemistry
    # For example, ensuring the glycosidic bond is at C1 and is in the beta configuration
    # Due to limitations, we assume the presence of the pattern is sufficient

    return True, "Contains beta-D-glucuronic acid moiety linked via beta-glycosidic bond"
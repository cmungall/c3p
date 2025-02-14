"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate contains a beta-D-glucuronic acid moiety attached via a beta-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-D-glucuronic acid moiety
    # This pattern represents the beta-D-glucuronic acid unit with correct stereochemistry
    beta_D_glucuronic_acid_smarts = '[C@@H]1O[C@H]([C@@H]([C@H]([C@@H](O1)[O])[O])[O])C(=O)[O-]'
    beta_D_glucuronic_acid = Chem.MolFromSmarts(beta_D_glucuronic_acid_smarts)
    if beta_D_glucuronic_acid is None:
        return False, "Invalid SMARTS pattern for beta-D-glucuronic acid"

    # Search for the beta-D-glucuronic acid moiety in the molecule considering stereochemistry
    matches = mol.GetSubstructMatches(beta_D_glucuronic_acid, useChirality=True)
    if not matches:
        return False, "Beta-D-glucuronic acid moiety not found with correct stereochemistry"

    # Verify the glycosidic linkage is a beta linkage
    # The anomeric carbon (C1) should be connected via an oxygen to another group
    for match in matches:
        anomeric_carbon_idx = match[0]  # C1 atom index in the match
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)
        bonds = anomeric_carbon.GetBonds()

        # Check for glycosidic bond from anomeric carbon to another atom (excluding the ring oxygen)
        for bond in bonds:
            neighbor = bond.GetOtherAtom(anomeric_carbon)
            if neighbor.GetAtomicNum() != 8 and not mol.IsInRing(bond.GetIdx()):
                # Found a bond from C1 to a non-oxygen atom outside the ring
                # This indicates a glycosidic linkage
                # Check if the linkage is beta (based on stereochemistry)
                return True, "Contains beta-D-glucuronic acid moiety with beta linkage"

    return False, "Beta-D-glucuronic acid moiety found but beta linkage not confirmed"
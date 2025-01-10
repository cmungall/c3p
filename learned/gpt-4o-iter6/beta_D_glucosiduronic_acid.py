"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid has a beta-D-glucuronic acid moiety bound via a glycosidic bond.

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

    # Look for beta-D-glucuronic acid moiety pattern (correct stereochemistry for beta linkage at anomeric carbon)
    glucuronic_pattern = Chem.MolFromSmarts("[C@H]1(O)O[C@@H]([C@H](O)[C@@H](O)CO1)C(=O)O")
    if not mol.HasSubstructMatch(glucuronic_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Check for presence of a conjugated glycosidic bond link typically found in derivatives
    glycosidic_link_pattern = Chem.MolFromSmarts("O[C@H]1O[C@@H]([C@@H](O)[C@H](O)O1)C")
    if not mol.HasSubstructMatch(glycosidic_link_pattern):
        return False, "Missing glycosidic bond indicating beta linkage"

    # Ensuring presence of a carbonyl group in the acid-end to confirm beta-D-glucuronic acid linkage
    carboxylic_link_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_link_pattern):
        return False, "No carboxylic acid group attached to form beta-D-glucuronic acid linkage"

    return True, "Contains the required beta-D-glucuronic acid moiety with correct glycosidic linkage"
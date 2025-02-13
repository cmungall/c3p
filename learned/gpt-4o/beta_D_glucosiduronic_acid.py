"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid involves a beta-D-glucuronic acid glycosidically linked to another moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a molecule is beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Pattern for glucuronic acid with appropriate beta linkage
    glucuronic_acid_pattern = Chem.MolFromSmarts("OC1[C@@H](O)[C@H](O)[C@H](O)[C@@H](O1)C(=O)O")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid substructure found"
    
    # Pattern for glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("O-C")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bond_matches:
        return False, "No glycosidic bond found"

    return True, "Contains beta-D-glucuronic acid with glycosidic linkage"

# Examples
# SMILES examples provided in the task description can be classified using the above function
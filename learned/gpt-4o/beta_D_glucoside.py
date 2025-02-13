"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is identified by the presence of a beta-D-glucose moiety linked through an ether bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-D-glucoside, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved SMARTS pattern for a beta-D-glucoside structure
    # Note: The stereochemistry in SMARTS might need further refinement depending on actual test results.

    # This pattern aims to capture: - a 6-membered pyranose ring
    #                               - hydroxyl group positions that are characteristic of glucose
    #                               - beta anomeric carbon (-O linkage)
    beta_d_glucoside_patterns = [
        Chem.MolFromSmarts("[C@@H]1(O[C@@H](CO)O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O)"),
        Chem.MolFromSmarts("O[C@@H]1[C@@H](O[C@H](O)[C@H](O)[C@H]1O)CO"),
        Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O1)CO")
    ]
    
    # Check if the molecule has the beta-D-glucoside pattern
    for pattern in beta_d_glucoside_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-D-glucoside substructure"
    
    return False, "No beta-D-glucoside substructure found"

# Example of usage:
# result, reason = is_beta_D_glucoside("OC[C@H]1O[C@@H](Oc2ccc3ccc(c2)-c2ccccc3c2)C1O")
# print(result, reason)
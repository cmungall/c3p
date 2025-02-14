"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for D-glucopyranosyl unit
    d_glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O")
    
    # Check if the molecule contains the D-glucose pattern
    if not mol.HasSubstructMatch(d_glucose_pattern):
        return False, "No D-glucose-derived glycoside group found"
    
    # Check for glycosidic linkage (presence of O-glycoside bond)
    # This might be more complex based on specific linkage types
    glycoside_pattern = Chem.MolFromSmarts("[C,O]OC[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycoside linkage found"

    return True, "Contains D-glucose-derived glycoside group with appropriate glycosidic linkage"

# Example usage:
# smiles = "OCC(CO)O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"  # 2-O-(alpha-D-glucopyranosyl)glycerol
# is_D_glucoside(smiles)
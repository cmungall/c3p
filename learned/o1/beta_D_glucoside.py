"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is any D-glucoside in which the anomeric centre has beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-D-glucoside
    # This pattern represents a glucose ring with beta linkage at C1
    beta_d_glucoside_smarts = "[C@H]1([O][#6])[O][C@@H]([C@H]([C@H]([C@H]1O)O)O)O"
    beta_d_glucoside_pattern = Chem.MolFromSmarts(beta_d_glucoside_smarts)
    if beta_d_glucoside_pattern is None:
        return False, "Error in SMARTS pattern"

    # Search for matches
    matches = mol.GetSubstructMatches(beta_d_glucoside_pattern)
    if matches:
        return True, "Contains beta-D-glucoside unit"
    else:
        return False, "No beta-D-glucoside unit found"
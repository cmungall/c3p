"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is a D-galactose residue with beta-configuration at the anomeric carbon
    connected via a glycosidic bond to any aglycone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-D-galactoside with correct stereochemistry
    beta_D_galactoside_smarts = '[O][C@H]1[C@@H](O[!#1])[C@@H](O)[C@H](O)[C@H](O)[C@H]1O'
    beta_D_galactoside_mol = Chem.MolFromSmarts(beta_D_galactoside_smarts)
    if beta_D_galactoside_mol is None:
        return False, "Error in SMARTS pattern for beta-D-galactoside"

    # Check if the molecule contains the beta-D-galactoside substructure
    if mol.HasSubstructMatch(beta_D_galactoside_mol):
        return True, "Contains beta-D-galactoside substructure"
    else:
        return False, "Does not contain beta-D-galactoside substructure"
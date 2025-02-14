"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A glucosiduronate is characterized by a beta-D-glucuronic acid moiety deprotonated
    at the carboxy group typically connected via an O-linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for beta-D-glucuronic acid-based structure
    # Look for cyclic structure with multiple hydroxyls and deprotonated carboxylate
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([O-]1)O)O)O)O)C(=O)[O-]")
    
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"
    
    # Ensure presence of O-linkage from the glucuronic moiety to another heteroatom
    o_linkage_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]")
    if not mol.HasSubstructMatch(o_linkage_pattern):
        return False, "No beta-D-glucuronic O-linkage found"

    # If all checks pass, it is a beta-D-glucosiduronate
    return True, "Contains beta-D-glucuronic acid moiety with proper O-linkage and deprotonation"
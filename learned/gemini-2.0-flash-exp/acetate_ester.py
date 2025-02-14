"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is a carboxylic ester where the carboxylic acid component is acetic acid.
    That means that the substructure must be CH3C(=O)O- attached to any other atom via the O
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for acetate ester group (CH3C(=O)O-)
    # The O atom is connected to any other atom.
    acetate_pattern = Chem.MolFromSmarts("CC(=O)O*")

    # Search for the substructure
    matches = mol.GetSubstructMatches(acetate_pattern)

    # Check if at least one match is found
    if len(matches) > 0:
        return True, "Contains an acetate group"
    else:
        return False, "Does not contain an acetate group"
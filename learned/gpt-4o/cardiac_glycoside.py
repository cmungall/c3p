"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import Lipinski

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the typical steroid backbone with four fused rings
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4=CC(CC4C3CCC21)C')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define a lactone ring attached to the steroid backbone
    lactone_pattern = Chem.MolFromSmarts('C1=CC(=O)OC1')
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for sugar moieties (glycosidic bonds)
    sugar_pattern = Chem.MolFromSmarts('[OH]C(C[OH])O')
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moieties found"

    # Check that the molecule has at least one sugar attachment to a steroid scaffold
    if not Lipinski.NumHDonors(mol) >= 3:
        return False, "Insufficient potential sugar attachments"

    return True, "Contains steroid backbone with lactone ring and glycosidic bond characteristic of cardiac glycosides"
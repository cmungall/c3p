"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a molecule containing a D-ribose sugar and a nucleobase.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the D-ribose SMARTS pattern with specific stereochemistry.
    #  C1(O)[C@H](O)[C@@H](CO)[C@H](O)O1
    ribose_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@@H]([O])[C@@H](CO)[C@H]([O])O1")
    if not mol.HasSubstructMatch(ribose_pattern):
         return False, "Molecule does not contain the correct D-ribose sugar"
    
    # check for nitrogen attachment to ribose carbon, meaning nucleobase attached. 
    base_attach_pattern = Chem.MolFromSmarts("[NX3][CX4]1[C@H]([O])[C@@H]([O])[C@@H](CO)[C@H]([O])O1")
    if not mol.HasSubstructMatch(base_attach_pattern):
        return False, "Molecule does not contain a nucleobase attached to the ribose sugar"

    return True, "Molecule contains a D-ribose sugar and a nucleobase attached to it."
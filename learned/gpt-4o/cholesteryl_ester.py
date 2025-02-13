"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined as a sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid nucleus pattern (generic for cholesterol backbone)
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
        
    # Define ester linkage pattern
    ester_linkage_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H0])]") # Pattern for ester C(=O)O
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found"

    # Optionally check for a fatty acid aliphatic chain (long carbon chain after ester linkage)
    # Although recognizing specific fatty acids is complex and may require inspecting specific atom counts and bond arrangements.
    
    return True, "Contains steroid backbone of cholesterol and ester linkage indicative of cholesteryl ester"
"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    A 3-substituted propionyl-CoA(4-) is an acyl-CoA(4-) oxoanion with a 3-substituted propionyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA backbone pattern
    # CoA backbone includes a pantothenic acid moiety, a phosphate group, and a diphosphate group
    coa_backbone_pattern = Chem.MolFromSmarts("[CX4][CX4]([OX2])([CX4])[CX4](=[OX1])[NX3][CX4][CX4](=[OX1])[NX3][CX4][CX4](=[OX1])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[PX4](=[OX1])([OX1-])[OX2][PX4](=[OX1])([OX1-])[OX2][CX4][CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])([CX4])[CX4]([OX2])(
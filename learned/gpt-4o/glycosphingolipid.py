"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid must contain both a sphingoid/ceramide backbone and a carbohydrate residue
    attached via a glycosidic linkage to O-1 of the sphingoid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General glycosidic linkage (CHO linkages)
    glycosidic_pattern = Chem.MolFromSmarts("[C][O][*R]") # C linked to an O part of a ring
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # General sphingoid base (e.g., base structures of ceramide)
    sphingoid_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@@H](O)CCCCCCCCCCC")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid/ceramide backbone found"

    # Both the patterns should be present for the classification as glycosphingolipid
    return True, "Contains sphingoid/ceramide backbone with glycosidic linkage to carbohydrate"
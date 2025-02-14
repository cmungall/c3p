"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: CHEBI:17704 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as 'Any member of class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some substances classified as bacterial glycolipids have the sugar part acylated by one or more fatty acids and the glycerol part may be absent.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (-O-C-C-O-) or fatty acid chains
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]")
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_fatty_acids = len(mol.GetSubstructMatches(fatty_acid_pattern)) >= 2
    
    # Look for saccharide (ring of carbons with multiple oxygens)
    saccharide_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")
    has_saccharide = mol.HasSubstructMatch(saccharide_pattern)
    
    # Check for glycosidic linkage (O-C-O-C)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2][CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    # Classify as glycolipid if it has a glycerol/fatty acid part linked to a saccharide via glycosidic bond
    if (has_glycerol or has_fatty_acids) and has_saccharide and glycosidic_matches:
        return True, "Contains glycerol/fatty acid part linked to a saccharide via glycosidic bond"
    else:
        return False, "Does not match the definition of a glycolipid"
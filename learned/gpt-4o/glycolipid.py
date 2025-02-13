"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as a 1,2-di-O-acylglycerol with a carbohydrate
    part joined via a glycosidic linkage.

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
    
    # Check for glycerol backbone with two acyl chains
    glycerol_pattern = Chem.MolFromSmarts("[C](O[C]=O)(O[C]=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 1,2-diacylglycerol structure detected"
    
    # Check for glycosidic linkage (sugar motif)
    # Monosaccharide pattern as Cn(O)Ox multiple hydroxyls
    sugar_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No glycosidic linkage to a sugar moiety found"
    
    # Count fatty acid chains (long carbon chains attached via esters or amides)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCC(=O)")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Expected at least 2 fatty acid chains, found {len(fatty_acid_matches)}"
    
    return True, "Structure matches a glycolipid with a glycerol backbone, acyl chains, and glycosidic linkage"
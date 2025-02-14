"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:87693 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is a fatty acid ester obtained by condensation of the carboxy group
    of tetradecanoic acid (myristic acid) with a hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for tetradecanoic acid fragment
    tetradecanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    tetradecanoic_acid_matches = mol.GetSubstructMatches(tetradecanoic_acid_pattern)
    
    # Look for ester bond (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check if tetradecanoic acid is esterified
    tetradecanoic_acid_esterified = False
    for match in tetradecanoic_acid_matches:
        acid_oxygen = match[-2]
        for ester_match in ester_matches:
            if acid_oxygen in ester_match:
                tetradecanoic_acid_esterified = True
                break
        if tetradecanoic_acid_esterified:
            break
            
    if not tetradecanoic_acid_esterified:
        return False, "Tetradecanoic acid not esterified"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 16:
        return False, "Too few carbons for tetradecanoate ester"
    if o_count < 3:
        return False, "Too few oxygens for tetradecanoate ester"
    
    return True, "Contains tetradecanoic acid esterified to an alcohol or phenol"
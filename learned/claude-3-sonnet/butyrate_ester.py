"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is any carboxylic ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Look for butyrate ester pattern: -O-C(=O)-CH2-CH2-CH3
    # This matches the butyric acid component connected via an ester bond
    butyrate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CH2X4][CH2X4][CH3X4]")
    matches = mol.GetSubstructMatches(butyrate_pattern)
    
    if not matches:
        return False, "No butyrate ester group found"
    
    # For each potential match, verify it's a butyrate ester by checking:
    # 1. The carbon chain is exactly 3 carbons (plus the carbonyl carbon)
    # 2. No branching in the butyric acid part
    for match in matches:
        # Get the atoms involved in the match
        ester_o, carbonyl_c, c1, c2, c3 = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Check that the carbons in the chain have correct number of hydrogens
        # C1 and C2 should be CH2, C3 should be CH3
        if (c1.GetTotalNumHs() == 2 and 
            c2.GetTotalNumHs() == 2 and 
            c3.GetTotalNumHs() == 3):
            return True, "Contains butyrate ester group (-O-C(=O)-CH2-CH2-CH3)"
            
    return False, "Found similar structure but not a true butyrate ester"
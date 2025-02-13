"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:17977 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (three-carbon chain) with two acyl, alkyl, or alk-1-enyl linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for a flexible glycerol backbone
    # Allowing for any connectivity around the central middle carbon
    glycerol_patterns = [
        Chem.MolFromSmarts("[CH2][CH]([OH])[CH2O]"),
        Chem.MolFromSmarts("[CH2][CH]([OH])[CH](O)*"),   # Allow for another hydroxyl or substituent
        Chem.MolFromSmarts("[CH]([OH])[CH2O][CH2]")      # Different arrangement
    ]

    # Check for a glycerol backbone using defined patterns
    glycerol_found = any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns)
    if not glycerol_found:
        return False, "No glycerol backbone found"

    # Look for the count of acyl, alkyl, alk-1-enyl groups (esters or ethers)
    # Any two variety of linkages at two positions - implementing flexible linkage pattern
    linkage_patterns = [
        Chem.MolFromSmarts("[O][CX3](=[O])"),  # ester linkage
        Chem.MolFromSmarts("[O][CH2][C]"),     # ether/alkyl linkage
        Chem.MolFromSmarts("[CH]=[CH]-[C]"),   # alk-1-enyl pattern
    ]

    linkage_count = sum(len(mol.GetSubstructMatches(pattern)) for pattern in linkage_patterns)
    if linkage_count < 2:
        return False, f"Found only {linkage_count} alkyl/acyl/alk-1-enyl linkages, need at least 2"

    return True, "Contains glycerol backbone with two acyl, alkyl, or alk-1-enyl linkages"